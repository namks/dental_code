################################################################################
# R code to perform TWAS
#
# Steps:
# 1. Run Main TWAS
# 2. Filter Significant Results & Annotate
################################################################################

suppressMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(biomaRt)
  library(parallel)
  library(doParallel)
  library(R.utils)
})

# --- Configuration ---
cfg <- list(
  sumstats_dir = "./SUMSTATS/",
  fusion_twas_script = "./fusion_twas-master/FUSION.assoc_test.R",
  weights_dir = "./fusion_twas-master/WEIGHTS/",
  ld_ref_prefix = "./fusion_twas-master/LDREF/1000G.EUR.",
  twas_output_dir = "./results/",
  num_cores = 10,
  chromosomes = 1:22,
  weights_pattern = "^GTExv8\\.EUR\\..*nofilter.*\\.pos$"
)

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
  flush.console()
}

# ---------------------------------------------------------------------------- #
# Step 1: Main TWAS Analysis
# ---------------------------------------------------------------------------- #
run_step_3_main_twas <- function(cfg) {
  log_message("--- [Step 1: Main TWAS Analysis] ---")
  
  sumstat_files <- list.files(cfg$sumstats_dir, pattern = "_chr[0-9]+\\.sumstats$", full.names = FALSE)
  disease_list <- unique(gsub("_chr[0-9]+\\.sumstats$", "", sumstat_files))
  
  weights_files <- list.files(cfg$weights_dir, pattern = cfg$weights_pattern, full.names = FALSE)
  tissue_list <- gsub("GTExv8\\.EUR\\.|\\.nofilter\\.pos$", "", weights_files)
  
  cl <- makeCluster(cfg$num_cores, outfile = "")
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("cfg", "disease_list", "log_message"), envir = environment())
  
  foreach(tissue = tissue_list, .packages = c("data.table", "R.utils")) %dopar% {
    for (disease in disease_list) {
      output_dir <- file.path(cfg$twas_output_dir, disease)
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      final_output_file <- file.path(output_dir, paste0(disease, "_", tissue, "_TWAS_combined_results_nofilter.txt"))
      if (file.exists(final_output_file)) next
      
      # Run per chromosome
      for (chr in cfg$chromosomes) {
        sumstats_file <- file.path(cfg$sumstats_dir, paste0(disease, "_chr", chr, ".sumstats"))
        weights_file <- file.path(cfg$weights_dir, paste0("GTExv8.EUR.", tissue, ".nofilter.pos"))
        output_file <- file.path(output_dir, paste0(disease, ".", tissue, ".chr", chr, ".dat"))
        
        if (!file.exists(sumstats_file)) next
        if (file.exists(output_file)) next
        
        cmd <- paste(
          "Rscript", shQuote(cfg$fusion_twas_script),
          "--sumstats", shQuote(sumstats_file),
          "--weights", shQuote(weights_file),
          "--weights_dir", shQuote(cfg$weights_dir),
          "--ref_ld_chr", shQuote(cfg$ld_ref_prefix),
          "--chr", chr,
          "--out", shQuote(output_file)
        )
        system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
      }
      
      # Merge Chromosomes
      chr_outputs <- list()
      for (chr in cfg$chromosomes) {
        base <- paste0(disease, ".", tissue, ".chr", chr)
        files <- c(file.path(output_dir, paste0(base, ".dat")), 
                   file.path(output_dir, paste0(base, ".dat.MHC")))
        for (f in files[file.exists(files)]) {
          dt <- fread(f)
          if (nrow(dt) > 0) chr_outputs[[length(chr_outputs) + 1]] <- dt
        }
      }
      if (length(chr_outputs) > 0) {
        fwrite(rbindlist(chr_outputs, fill = TRUE), final_output_file, sep = "\t")
      }
    }
  }
  stopCluster(cl)
  log_message("--- [Step 1 Complete] ---")
}

# ---------------------------------------------------------------------------- #
# Step 2: Filter & Annotate
# ---------------------------------------------------------------------------- #
run_step_4_filter_annotate <- function(cfg) {
  log_message("--- [Step 2: Filter TWAS Results] ---")
  
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  disease_dirs <- list.dirs(cfg$twas_output_dir, recursive = FALSE, full.names = TRUE)
  disease_dirs <- disease_dirs[!basename(disease_dirs) %in% c("PERMUTATION", "COLOC", "top_files")]
  
  for (disease_dir in disease_dirs) {
    disease_name <- basename(disease_dir)
    files <- list.files(disease_dir, pattern = "nofilter.*\\.txt$", full.names = TRUE)
    results <- list()
    
    for (f in files) {
      tryCatch({
        df <- fread(f)
        tissue_name <- basename(f) %>% sub(paste0(disease_name, "_"), "", .) %>% sub("_TWAS_combined_results.*", "", .)
        
        # Bonferroni Filter
        threshold <- 0.05 / max(1, (nrow(df) - 1))
        df_sig <- df[TWAS.P < threshold]
        
        if (nrow(df_sig) > 0) {
          df_sig[, Tissue := tissue_name]
          results[[length(results) + 1]] <- df_sig
        }
      }, error = function(e) return(NULL))
    }
    
    if (length(results) > 0) {
      combined <- rbindlist(results, fill=TRUE)
      
      # Annotation
      ids_clean <- unique(gsub("\\..*", "", combined$ID))
      gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ids_clean, mart = mart)
      
      combined[, ID_clean := gsub("\\..*", "", ID)]
      final <- merge(combined, gene_map, by.x="ID_clean", by.y="ensembl_gene_id", all.x=TRUE)
      final[, ID_clean := NULL]
      
      out_file <- file.path(cfg$twas_output_dir, paste0(disease_name, "_All_Tissues_TWAS_Filtered_Results.txt"))
      fwrite(final, out_file, sep = "\t")
      log_message(paste("Saved Filtered Results:", out_file))
    }
  }
  log_message("--- [Step 2 Complete] ---")
}

# Execution
if (!interactive()) {
  run_step_3_main_twas(cfg)
  run_step_4_filter_annotate(cfg)
}