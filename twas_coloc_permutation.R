################################################################################
# R code to perform COLOC & permutation test
#
# Steps:
# 6. Run Permutation Test
# 7. Run COLOC Analysis - use coloc 3.2.1
# 8. Merge COLOC Results
################################################################################

suppressMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(biomaRt)
  library(parallel)
  library(doParallel)
})

cfg <- list(
  sumstats_split_dir = "./SUMSTATS/split/",
  fusion_twas_script = "./fusion_twas-master/FUSION.assoc_test.R",
  weights_dir = "./fusion_twas-master/WEIGHTS/",
  ld_ref_prefix = "./fusion_twas-master/LDREF/1000G.EUR.",
  panel_n_file = "./GTEx_sample_sizes.txt", # contains PANEL, N
  gwas_n_file = "./gwas_sample_sizes.txt",  # contains Disease, GWAS_N
  
  # Input: step 02 Output
  joint_results_file = "./results/top_files/all_joint_included_results_for_post.txt",
  
  perm_output_dir = "./results/PERMUTATION/",
  coloc_output_dir = "./results/COLOC/",
  coloc_merged_output_file = "./results/COLOC/merged_coloc_results_annotated.txt",
  
  num_cores = 100,
  perm_n = 1000
)

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
  flush.console()
}

# ---------------------------------------------------------------------------- #
# Step 6: Permutation test and COLOC Analysis
# ---------------------------------------------------------------------------- #
run_step_5_perm_coloc <- function(cfg) {
  log_message("--- [Step 5: Running Permutation & COLOC on Joint Results] ---")
  
  if (!file.exists(cfg$joint_results_file)) stop("Joint results file not found. Run Post-process first.")
  if (!file.exists(cfg$gwas_n_file)) stop("GWAS sample size file not found.")
  
  joint_df <- fread(cfg$joint_results_file)
  gwas_n_df <- fread(cfg$gwas_n_file)
  
  if (nrow(joint_df) == 0) {
    log_message("No significant joint results to process.")
    return()
  }
  
  # Cluster Setup
  cl <- makeCluster(cfg$num_cores, outfile = "")
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("cfg", "joint_df", "gwas_n_df", "log_message"), envir = environment())
  
  foreach(i = 1:nrow(joint_df), .packages = c("data.table", "stringr")) %dopar% {
    row <- joint_df[i, ]
    
    # Extract Info
    gene_id <- row$ID
    tissue <- row$Tissue      
    disease_full <- row$Disease 
    chr <- row$CHR
    
    # Get GWAS N
    split_name <- strsplit(disease_full, "_merged")[[1]][1]
    gwas_key <- gsub("step2_", "", split_name)
    GWAS_N <- gwas_n_df[Disease == gwas_key, GWAS_N]
    
    if (length(GWAS_N) == 0) return(NULL)
    
    sumstats_file <- file.path(cfg$sumstats_split_dir, paste0(disease_full, "_chr", chr, ".sumstats"))
    weights_file <- file.path(cfg$weights_dir, paste0("GTExv8.EUR.", tissue, ".nofilter.pos"))
    
    if (!file.exists(sumstats_file) || !file.exists(weights_file)) return(NULL)
    
    # --- Permutation ---
    perm_dir <- file.path(cfg$perm_output_dir, disease_full)
    if (!dir.exists(perm_dir)) dir.create(perm_dir, recursive = TRUE)
    
    perm_file <- file.path(perm_dir, paste0(disease_full, ".", tissue, ".chr", chr, ".", gene_id, ".perm"))
    
    if (!file.exists(perm_file)) {
      cmd <- paste(
        "Rscript", shQuote(cfg$fusion_twas_script),
        "--sumstats", shQuote(sumstats_file),
        "--weights", shQuote(weights_file),
        "--weights_dir", shQuote(cfg$weights_dir),
        "--ref_ld_chr", shQuote(cfg$ld_ref_prefix),
        "--chr", chr,
        "--out", shQuote(perm_file),
        "--perm", cfg$perm_n,
        "--PANELN", shQuote(cfg$panel_n_file),
        "--GWASN", GWAS_N
      )
      system(cmd, ignore.stdout = TRUE)
    }
    
    # --- COLOC ---
    coloc_dir <- file.path(cfg$coloc_output_dir, disease_full)
    if (!dir.exists(coloc_dir)) dir.create(coloc_dir, recursive = TRUE)
    
    coloc_file <- file.path(coloc_dir, paste0(disease_full, ".", tissue, ".chr", chr, ".", gene_id, ".coloc.dat"))
    
    if (!file.exists(coloc_file)) {
      cmd <- paste(
        "Rscript", shQuote(cfg$fusion_twas_script),
        "--sumstats", shQuote(sumstats_file),
        "--weights", shQuote(weights_file),
        "--weights_dir", shQuote(cfg$weights_dir),
        "--ref_ld_chr", shQuote(cfg$ld_ref_prefix),
        "--chr", chr,
        "--out", shQuote(coloc_file),
        "--coloc_P", 1,
        "--PANELN", shQuote(cfg$panel_n_file),
        "--GWASN", GWAS_N
      )
      system(cmd, ignore.stdout = TRUE)
    }
  }
  
  stopCluster(cl)
  log_message("--- [Step 7 Complete] ---")
}

# ---------------------------------------------------------------------------- #
# Step 8: Merge COLOC Results
# ---------------------------------------------------------------------------- #
run_step_8_merge_coloc <- function(cfg) {
  log_message("--- [Step 7: Scanning & Merging All Existing Results] ---")
  
  files <- list.files(cfg$coloc_output_dir, pattern = "\\.coloc\\.dat(\\.MHC)?$", recursive = TRUE, full.names = TRUE)
  
  if (length(files) == 0) {
    log_message("No COLOC files found in directory.")
    return()
  }
  
  log_message(sprintf("Found %d COLOC/MHC files. Starting processing...", length(files)))
  
  res_list <- list()
  
  for (f in files) {
    tryCatch({
      fname <- basename(f)
      disease_name <- basename(dirname(f)) 
      
      ensg <- str_extract(fname, "ENSG[0-9]{11}(\\.\\d+)?")
      
      rem <- paste0(disease_name, ".")
      if (startsWith(fname, rem)) {
        rest <- sub(rem, "", fname, fixed = TRUE)
        parts <- strsplit(rest, "\\.")[[1]]
        chr_idx <- grep("^chr\\d+", parts)[1]
        tissue_val <- paste(parts[1:(chr_idx-1)], collapse = ".")
      } else {
        tissue_val <- unlist(strsplit(fname, "\\."))[2]
      }
      
      dt <- fread(f)
      if (nrow(dt) == 0) next
      
      # --- C. Filter ID (COLOC) ---
      dt_sub <- dt[ID == ensg]
      
      if (nrow(dt_sub) == 0) {
        ensg_clean <- str_replace(ensg, "\\.\\d+$", "")
        dt_sub <- dt[str_replace(ID, "\\.\\d+$", "") == ensg_clean]
      }
      
      if (nrow(dt_sub) > 0) {
        
        dt_sub <- dt_sub[1]
        
        dt_sub[, Disease := disease_name]
        dt_sub[, Tissue := tissue_val]
        dt_sub[, ID := ensg] 
        
        if (endsWith(fname, ".MHC")) {
          dt_sub[, Note := "MHC_File"]
        } else {
          dt_sub[, Note := "Standard"]
        }
        
        # --- D. Merge Permutation ---
        perm_fname <- sub("\\.coloc\\.dat(\\.MHC)?$", ".perm\\1", fname)
        perm_path <- file.path(cfg$perm_output_dir, disease_name, perm_fname)
        
        dt_sub[, `:=`(PERM.PV = as.numeric(NA), 
                      PERM.N = as.numeric(NA), 
                      PERM.ANL_PV = as.numeric(NA))]
        
        if (file.exists(perm_path)) {
          perm_dt <- fread(perm_path)
          
          if (nrow(perm_dt) > 0) {
            curr_id <- dt_sub$ID[1]
            
            target_row <- perm_dt[ID == curr_id]
            
            if (nrow(target_row) == 0) {
              curr_id_clean <- str_replace(curr_id, "\\.\\d+$", "")
              perm_dt[, ID_clean := str_replace(ID, "\\.\\d+$", "")]
              target_row <- perm_dt[ID_clean == curr_id_clean]
            }
            
            if (nrow(target_row) > 0) {
              if ("PERM.PV" %in% names(target_row)) dt_sub[, PERM.PV := as.numeric(target_row$PERM.PV[1])]
              if ("PERM.N" %in% names(target_row))  dt_sub[, PERM.N := as.numeric(target_row$PERM.N[1])]
              if ("PERM.ANL_PV" %in% names(target_row)) dt_sub[, PERM.ANL_PV := as.numeric(target_row$PERM.ANL_PV[1])]
            }
          }
        }
        
        res_list[[length(res_list) + 1]] <- dt_sub
      }
      
    }, error = function(e) { return(NULL) })
    
    if (length(res_list) %% 1000 == 0) cat(sprintf("...processed %d files\n", length(res_list)))
  }
  
  if (length(res_list) == 0) {
    log_message("No valid data extracted from scanned files.")
    return()
  }
  
  merged_coloc <- rbindlist(res_list, fill = TRUE)
  log_message(sprintf("Total merged rows: %d", nrow(merged_coloc)))
  
  if (file.exists(cfg$joint_results_file)) {
    joint_df <- fread(cfg$joint_results_file)
    joint_cols <- c("ID", "Tissue", "Disease", "JOINT.BETA", "JOINT.BETA.SE", "JOINT.Z", "JOINT.P")
    use_cols <- intersect(names(joint_df), joint_cols)
    
    merged_coloc[, `:=`(ID=as.character(ID), Tissue=as.character(Tissue), Disease=as.character(Disease))]
    joint_subset <- joint_df[, ..use_cols]
    joint_subset[, `:=`(ID=as.character(ID), Tissue=as.character(Tissue), Disease=as.character(Disease))]
    
    merged_coloc <- merge(merged_coloc, joint_subset, by = c("ID", "Tissue", "Disease"), all.x = TRUE)
  }
  
  # --- ANNOTATION Using biomaRt ---
  log_message("Running biomaRt annotation...")
  
  tryCatch({
    merged_coloc[, ID_clean := str_replace(ID, "\\.\\d+$", "")]
    unique_ids <- unique(merged_coloc$ID_clean)
    
    if(length(unique_ids) > 0) {
      mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      
      gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        filters = "ensembl_gene_id",
                        values = unique_ids,
                        mart = mart)
      
      setDT(gene_map)
      
      if ("hgnc_symbol" %in% names(merged_coloc)) merged_coloc[, hgnc_symbol := NULL]
      
      merged_coloc <- merge(merged_coloc, gene_map, 
                            by.x = "ID_clean", by.y = "ensembl_gene_id", 
                            all.x = TRUE)
      
      log_message("biomaRt annotation completed.")
    } else {
      log_message("No IDs found for biomaRt annotation.")
    }
    
    merged_coloc[, ID_clean := NULL]
    
  }, error = function(e) {
    log_message(paste("Error during biomaRt annotation:", e$message))
    log_message("Proceeding without HGNC symbols.")
  })
  
  
  preferred_order <- c("Disease", "Tissue", "ID", "hgnc_symbol", "CHR", "P0", "P1", "P2", "P3", "P4", "PP.H4.abf", 
                       "PERM.PV", "PERM.N", "PERM.ANL_PV", 
                       "JOINT.BETA", "JOINT.BETA.SE", "JOINT.Z", "JOINT.P", "Note")
  
  existing_cols <- names(merged_coloc)
  col_order <- c(intersect(preferred_order, existing_cols), setdiff(existing_cols, preferred_order))
  setcolorder(merged_coloc, col_order)
  
  fwrite(merged_coloc, cfg$coloc_merged_output_file, sep = "\t", na = "NA", quote = FALSE)
  log_message(paste("✅ Final Merged Results saved to:", cfg$coloc_merged_output_file))
}

# Execution
if (!interactive()) {
  run_step_5_perm_coloc(cfg)
  run_step_8_merge_coloc(cfg)
}

