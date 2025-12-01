# R code to perform Integrated eQTL Mendelian Randomization (MR) Pipeline

# =============================================================================
# STEP 0: SETUP - LOAD LIBRARIES AND DEFINE PATHS
# =============================================================================

message("STEP 0: Loading required packages...")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(future)       
  library(future.apply) 
})

data.table::setDTthreads(1)
options(future.stdout = TRUE)
base_dir <- getwd()
plink_path <- file.path(base_dir, "plink", "plink.exe") 

# Input Data Paths
gtex_raw_dir     <- file.path(base_dir, "data", "gtex_v8_raw")
gtex_lookup_dir  <- file.path(base_dir, "data", "gtex_v8_lookup")
ld_ref_dir       <- file.path(base_dir, "data", "ld_reference")
gwas_data_dir    <- file.path(base_dir, "data", "gwas_summary_stats")

# Configuration Files
config_dir       <- file.path(base_dir, "config")
master_config_file <- file.path(config_dir, "merged_coloc_results_annotated.txt") #fill in the hgnc_symbol with custom_label for analysis.

# Output Results Paths
parsed_output_dir <- file.path(base_dir, "results", "1_gtex_parsed")
goi_output_dir    <- file.path(base_dir, "results", "2_goi_eqtl")
mr_output_dir     <- file.path(base_dir, "results", "3_mr_analysis")

# --- Create Output Directories ---
dir.create(parsed_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(goi_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mr_output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper Function ---
clean_gene_folder <- function(symbol) {
  symbol %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_$", "")
}

# --- Load and Filter Master Config File ---
message("STEP 0: Loading and filtering master configuration file...")
master_config_data <- fread(master_config_file)


# --- SETUP PARALLEL PROCESSING ---
plan(multisession, workers = 8) 
message(paste("STEP 0: Setup complete. Parallel processing enabled with 8 workers."))


# =============================================================================
# STEP 1: PARSE RAW GTEx v8 eQTL FILES
# =============================================================================
message("\n--- STEP 1: Parsing raw GTEx v8 eQTL files ---")
eqtl_files_raw <- list.files(gtex_raw_dir, pattern = "signif_variant_gene_pairs.txt$", full.names = TRUE)

invisible(future_lapply(eqtl_files_raw, function(eqtl_file) {
  tryCatch({
    tissue_name <- gsub("^.*[\\\\/]|\\.v8\\.signif_variant_gene_pairs\\.txt$", "", eqtl_file) 
    output_file <- file.path(parsed_output_dir, paste0(tissue_name, "_parsed.txt"))
    
    if (file.exists(output_file)) {
      message(paste("  Skipping:", basename(output_file), "(already exists)"))
      return(NULL) 
    }
    
    message(paste("  Parsing:", basename(eqtl_file)))
    eqtl_data <- fread(eqtl_file)
    
    var_split <- str_split_fixed(eqtl_data$variant_id, "_", 5)
    eqtl_data <- eqtl_data %>%
      mutate(
        chr = gsub("chr", "", var_split[,1]),
        pos = as.integer(var_split[,2]),
        ref = var_split[,3],
        alt = var_split[,4],
        MarkerID = paste(chr, pos, ref, alt, sep = ":")
      )
    
    fwrite(eqtl_data, output_file, sep = "\t")
    message(paste("  Saved:", basename(output_file)))
  }, error = function(e) {
    message(paste("  ERROR parsing", basename(eqtl_file), ":", e$message))
  })
}))
message("--- STEP 1: Parsing complete. ---")


# =============================================================================
# STEP 2: EXTRACT GENES OF INTEREST (GOI) FROM PARSED FILES
# =============================================================================
message("\n--- STEP 2: Extracting Genes of Interest (GOI) ---")
goi_map <- master_config_data %>%
  distinct(Disease, hgnc_symbol, ID, CHR) %>%
  rename(gene_symbol = hgnc_symbol, gene_id = ID, chr = CHR) %>%
  mutate(chr = as.character(chr))

parsed_files <- list.files(parsed_output_dir, pattern = "_parsed.txt$", full.names = TRUE)

invisible(future_lapply(1:nrow(goi_map), function(i) {
  tryCatch({
    gene_symbol <- goi_map$gene_symbol[i]
    gene_id_full <- goi_map$gene_id[i]
    chr_value <- goi_map$chr[i]
    gene_folder <- clean_gene_folder(gene_symbol)
    
    message(paste("\n  Processing GOI:", gene_symbol, "(ID:", gene_id_full, ", chr", chr_value, ")"))
    
    output_dir_base <- file.path(goi_output_dir, gene_folder, "filtered")
    if (!dir.exists(output_dir_base)) dir.create(output_dir_base, recursive = TRUE)
    
    for (parsed_file in parsed_files) {
      tissue_name <- gsub("^.*[\\\\/]|_parsed.txt$", "", parsed_file)
      output_file <- file.path(output_dir_base, paste0(tissue_name, "_filtered.txt"))
      
      if (file.exists(output_file)) {
        next
      }
      
      eqtl_data <- fread(parsed_file) %>% mutate(chr = as.character(chr))
      
      filtered_data <- eqtl_data %>%
        filter(gene_id == gene_id_full & chr == chr_value)
      
      if (nrow(filtered_data) == 0) next
      
      final_data <- filtered_data %>%
        select(
          gene_id, variant_id, MarkerID, chr, pos, ref, alt,
          beta = slope, se = slope_se, pval = pval_nominal, maf,
          tss_distance, ma_samples, ma_count, pval_nominal_threshold,
          min_pval_nominal, pval_beta
        )
      
      fwrite(final_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      message(paste("    Saved:", gene_symbol, "-", basename(output_file)))
    }
  }, error = function(e) {
    message(paste("  ERROR processing gene index", i, ":", e$message))
  })
}))
message("--- STEP 2: GOI extraction complete. ---")


# =============================================================================
# STEP 3: MAP rsIDs TO GOI DATA
# =============================================================================
message("\n--- STEP 3: Mapping rsIDs to GOI data (in parallel) ---")
gene_dirs <- list.dirs(goi_output_dir, recursive = FALSE, full.names = TRUE)

invisible(future_lapply(gene_dirs, function(gene_path) {
  tryCatch({
    gene_name <- basename(gene_path)
    input_dir <- file.path(gene_path, "filtered")
    output_dir <- file.path(gene_path, "filtered_with_rsid")
    
    if (!dir.exists(input_dir)) return(NULL)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    file_list <- list.files(input_dir, pattern = "_filtered.txt$", full.names = TRUE)
    
    for (file in file_list) {
      output_file <- file.path(output_dir, basename(file))
      if (file.exists(output_file)) next
      
      dat <- fread(file)
      if (nrow(dat) == 0) next
      
      snp_chr <- unique(dat$chr)
      if (length(snp_chr) != 1) {
        message(paste("  ERROR: Ambiguous CHR in:", file))
        next
      }
      
      lookup_file <- file.path(gtex_lookup_dir, paste0("lookup_chr", snp_chr, ".txt"))
      if (!file.exists(lookup_file)) {
        message(paste("  ERROR: Lookup file not found:", lookup_file))
        next
      }
      lookup <- fread(lookup_file, select = c("variant_id", "rs_id_dbSNP151_GRCh38p7"))
      
      if (!"variant_id" %in% names(dat)) {
        if (all(c("chr", "pos", "ref", "alt") %in% names(dat))) {
          dat[, variant_id := paste0("chr", chr, "_", pos, "_", ref, "_", alt, "_b38")]
        } else {
          message(paste("  ERROR: Cannot create variant_id, missing columns in:", file))
          next
        }
      }
      
      merged <- merge(dat, lookup, by = "variant_id", all.x = TRUE)
      setnames(merged, "rs_id_dbSNP151_GRCh38p7", "rsid")
      
      cleaned <- merged[!is.na(rsid) & rsid != "."]
      
      if (nrow(cleaned) > 0) {
        setcolorder(cleaned, c("gene_id", "rsid", setdiff(names(cleaned), c("gene_id", "rsid"))))
        fwrite(cleaned, output_file, sep = "\t")
        message(paste("  Mapped rsIDs for:", gene_name, "-", basename(output_file)))
      }
    }
  }, error = function(e) {
    message(paste("  ERROR in rsID mapping for", basename(gene_path), ":", e$message))
  })
}))
message("--- STEP 3: rsID mapping complete. ---")


# =============================================================================
# STEP 4: LD CLUMPING (using PLINK)
# =============================================================================
message("\n--- STEP 4: Performing LD Clumping with PLINK ---")
gene_dirs_clump <- list.dirs(goi_output_dir, recursive = FALSE, full.names = TRUE)

invisible(future_lapply(gene_dirs_clump, function(gene_path) {
  tryCatch({
    gene_name <- basename(gene_path)
    eqtl_dir_clump <- file.path(gene_path, "filtered_with_rsid")
    output_dir_clump <- file.path(gene_path, "clumped", "byrska")
    
    if (!dir.exists(eqtl_dir_clump)) return(NULL)
    if (!dir.exists(output_dir_clump)) dir.create(output_dir_clump, recursive = TRUE)
    
    eqtl_files_clump <- list.files(eqtl_dir_clump, pattern = "filtered.txt$", full.names = TRUE)
    
    for (eqtl_file in eqtl_files_clump) {
      tissue_name <- gsub("_filtered.txt$", "", basename(eqtl_file))
      message(paste("  Clumping:", gene_name, "-", tissue_name))
      
      eqtl_data <- fread(eqtl_file)
      if (!"rsid" %in% colnames(eqtl_data) || nrow(eqtl_data) == 0) {
        message("    WARNING: No rsIDs to clump.")
        next
      }
      
      chr_val <- unique(eqtl_data$chr)
      if (length(chr_val) != 1) {
        message("    ERROR: Ambiguous CHR for clumping.")
        next
      }
      
      output_prefix <- file.path(output_dir_clump, paste0("clumped_", tissue_name, "_chr", chr_val))
      final_clumped_file <- paste0(output_prefix, "_final.csv")
      
      if (file.exists(final_clumped_file)) {
        message("    Skipping (already exists).")
        next
      }
      
      plink_ref_panel <- file.path(ld_ref_dir, paste0("chr", chr_val, "_hg38_no_dups"))
      if (!file.exists(paste0(plink_ref_panel, ".bed"))) {
        message(paste("    ERROR: LD Reference panel not found:", plink_ref_panel, ".bed"))
        next
      }
      
      plink_input_file <- tempfile(fileext = ".txt")
      plink_input <- eqtl_data %>% select(rsid, pval)
      setnames(plink_input, c("rsid", "pval"), c("SNP", "P"))
      fwrite(plink_input, plink_input_file, sep = "\t")

      cmd <- paste0(
        '"', plink_path, '"', # Quoted PLINK executable path
        " --bfile ", shQuote(plink_ref_panel, type = "cmd"),
        " --clump ", shQuote(plink_input_file, type = "cmd"),
        " --clump-r2 0.1 ",
        " --clump-field P ",
        " --clump-snp-field SNP ",
        " --out ", shQuote(output_prefix, type = "cmd")
      )

      system(cmd, intern = TRUE)
      
      clumped_snps_file <- paste0(output_prefix, ".clumped")
      if (file.exists(clumped_snps_file)) {
        clumped_snps <- fread(clumped_snps_file)
        fwrite(clumped_snps, final_clumped_file, row.names = FALSE)
        message(paste("    SUCCESS: Saved", basename(final_clumped_file)))
      } else {
        message("    ERROR: PLINK clumping failed to produce output.")
      }
      
      # Clean up temporary files
      unlink(c(plink_input_file, paste0(output_prefix, c(".log", ".nosex", ".clumped"))))
    }
  }, error = function(e) {
    message(paste("  ERROR in clumping for", basename(gene_path), ":", e$message))
  })
}))
message("--- STEP 4: LD Clumping complete. ---")


# =============================================================================
# STEP 5: CREATE EXPOSURE DATA
# =============================================================================
message("\n--- STEP 5: Creating final MR Exposure (IV) datasets ---")
gene_dirs_exp <- list.dirs(goi_output_dir, recursive = FALSE, full.names = TRUE)

invisible(future_lapply(gene_dirs_exp, function(gene_path) {
  tryCatch({
    gene_name <- basename(gene_path)
    eqtl_dir_exp <- file.path(gene_path, "filtered_with_rsid")
    clumped_dir_exp <- file.path(gene_path, "clumped", "byrska")
    output_dir_exp <- file.path(gene_path, "mr_input")
    
    if (!dir.exists(eqtl_dir_exp) || !dir.exists(clumped_dir_exp)) return(NULL)
    if (!dir.exists(output_dir_exp)) dir.create(output_dir_exp, recursive = TRUE)
    
    clumped_files <- list.files(clumped_dir_exp, pattern = "_final\\.csv$", full.names = TRUE)
    
    for (clumped_path in clumped_files) {
      base <- gsub("^clumped_|_chr\\d+_final\\.csv$", "", basename(clumped_path))
      output_file <- file.path(output_dir_exp, paste0("mr_exposure_", base, ".txt"))
      
      if (file.exists(output_file)) next
      
      eqtl_path <- file.path(eqtl_dir_exp, paste0(base, "_filtered.txt"))
      if (!file.exists(eqtl_path)) {
        message(paste("  ERROR: eQTL file not found:", gene_name, "-", base))
        next
      }
      
      eqtl_data <- fread(eqtl_path)
      clumped_data <- fread(clumped_path)
      
      filtered_eqtl <- eqtl_data %>% filter(pos %in% clumped_data$BP)
      
      if (nrow(filtered_eqtl) > 0) {
        fwrite(filtered_eqtl, output_file, sep = "\t")
        message(paste("  Saved MR Exposure:", gene_name, "-", basename(output_file)))
      }
    }
  }, error = function(e) {
    message(paste("  ERROR creating exposure data for", basename(gene_path), ":", e$message))
  })
}))
message("--- STEP 5: Exposure data creation complete. ---")


# =============================================================================
# STEP 6: EXTRACT OUTCOME DATA
# =============================================================================
message("\n--- STEP 6: Extracting matching Outcome data ---")
mr_combinations <- master_config_data %>%
  filter(!is.na(hgnc_symbol), !is.na(ID)) %>%
  rename(gene_symbol = hgnc_symbol, gene_id = ID) %>%
  distinct(gene_id, gene_symbol, Disease)

invisible(future_lapply(1:nrow(mr_combinations), function(i) {
  tryCatch({
    gene_symbol_label <- mr_combinations$gene_symbol[i]
    gene_safe <- clean_gene_folder(mr_combinations$gene_symbol[i])
    disease <- mr_combinations$Disease[i]
    
    message(paste("  Extracting Outcome for:", gene_symbol_label, "vs", disease))
    
    exposure_dir <- file.path(goi_output_dir, gene_safe, "mr_input")
    if (!dir.exists(exposure_dir)) return(NULL)
    
    iv_files <- list.files(exposure_dir, pattern = "\\.txt$", full.names = TRUE)
    
    if (length(iv_files) == 0) {
      message(paste("    No IV files found for", gene_symbol_label))
      return(NULL)
    }
    
    all_ivs_for_gene <- lapply(iv_files, fread) %>%
      bind_rows() %>%
      distinct(rsid, chr, ref, alt)
    
    if (nrow(all_ivs_for_gene) == 0) return(NULL)
    
    outcome_output_dir <- file.path(mr_output_dir, "temp_outcome_data", gene_safe, disease)
    dir.create(outcome_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    matched_results <- list()
    
    for (chr_num in unique(all_ivs_for_gene$chr)) {
      gwas_file <- file.path(gwas_data_dir, paste0("step2_", disease, "_chr", chr_num, "_EUR.txt"))
      if (!file.exists(gwas_file)) {
        message(paste("  ERROR: GWAS file not found for", disease, "chr", chr_num))
        next
      }
      
      gwas_data <- fread(gwas_file)
      iv_data_chr <- all_ivs_for_gene %>% filter(chr == chr_num)
      
      direct_idx <- with(iv_data_chr, paste(rsid, ref, alt))
      flip_idx   <- with(iv_data_chr, paste(rsid, alt, ref))
      gwas_idx   <- with(gwas_data, paste(MarkerID, Allele1, Allele2))
      
      matched <- gwas_data[gwas_idx %in% direct_idx | gwas_idx %in% flip_idx, ]
      
      matched <- matched %>%
        filter(AF_Allele2 >= 0.01, AF_Allele2 <= 0.99)
      
      if (nrow(matched) > 0) {
        matched_results[[as.character(chr_num)]] <- matched
      }
    }
    
    if (length(matched_results) > 0) {
      tissue_result <- bind_rows(matched_results)
      out_path <- file.path(outcome_output_dir, "matched_outcome_SNPs.txt")
      fwrite(tissue_result, out_path, sep = "\t")
      message(paste("    Saved matched outcome SNPs for:", gene_symbol_label, "vs", disease))
    } else {
      message(paste("    WARNING: No matching outcome SNPs found for:", gene_symbol_label, "vs", disease))
    }
  }, error = function(e) {
    message(paste("  ERROR extracting outcome for gene index", i, ":", e$message))
  })
}))
message("--- STEP 6: Outcome data extraction complete. ---")


# =============================================================================
# STEP 7: PERFORM TWO-SAMPLE MR ANALYSIS
# =============================================================================
message("\n--- STEP 7: Performing Two-Sample MR Analysis ---")
mr_run_list <- master_config_data %>%
  filter(!is.na(hgnc_symbol)) %>%
  rename(gene_symbol = hgnc_symbol, gene_id = ID) %>%
  distinct(gene_symbol, gene_id, Disease, Tissue)

fmt_pval <- function(p_vector) {
  p_formatted <- ifelse(is.na(p_vector) | !is.numeric(p_vector), 
                        NA_character_, 
                        format(round(p_vector, 4), nsmall = 4))
  
  p_formatted <- ifelse(!is.na(p_vector) & p_vector < 0.001, 
                        "<0.001", 
                        p_formatted)
  
  return(p_formatted)
}

all_results_list <- future_lapply(1:nrow(mr_run_list), function(i) {
  tryCatch({
    gene_symbol_label <- mr_run_list$gene_symbol[i]
    gene_safe <- clean_gene_folder(mr_run_list$gene_symbol[i])
    disease <- mr_run_list$Disease[i]
    tissue <- mr_run_list$Tissue[i]
    
    tissue_file_id <- str_replace_all(tissue, " ", "_")
    
    message(paste("  MR RUN:", gene_symbol_label, "(", gene_safe, ") -", tissue, "VS", disease))
    
    exposure_file <- file.path(goi_output_dir, gene_safe, "mr_input", paste0("mr_exposure_", tissue_file_id, ".txt"))
    outcome_file  <- file.path(mr_output_dir, "temp_outcome_data", gene_safe, disease, "matched_outcome_SNPs.txt")
    final_result_dir <- file.path(mr_output_dir, "final_results", gene_safe, "filtered")
    
    if (!file.exists(exposure_file)) { message("    ERROR: Exposure file not found."); return(NULL) }
    if (!file.exists(outcome_file)) { message("    ERROR: Outcome file not found."); return(NULL) }
    if (!dir.exists(final_result_dir)) dir.create(final_result_dir, recursive = TRUE)
    
    exp_raw <- fread(exposure_file)
    exp <- exp_raw %>%
      rename(SNP = rsid, effect_allele.exposure = alt, other_allele.exposure = ref,
             beta.exposure = beta, se.exposure = se) %>%
      mutate(
        eaf.exposure = if ("maf" %in% colnames(exp_raw)) exp_raw$maf else NA,
        pval.exposure = if ("pval" %in% colnames(exp_raw)) exp_raw$pval else NA,
        id.exposure = paste0(gene_safe, "_", tissue),
        exposure = paste0(gene_symbol_label, " expression in ", tissue)
      ) %>%
      select(all_of(c("SNP", "beta.exposure", "se.exposure", "eaf.exposure",
                      "effect_allele.exposure", "other_allele.exposure",
                      "pval.exposure", "id.exposure", "exposure")))
    
    out_raw <- fread(outcome_file)
    out <- out_raw %>%
      rename(SNP = MarkerID, effect_allele.outcome = Allele2, other_allele.outcome = Allele1,
             beta.outcome = BETA, se.outcome = SE, pval.outcome = p.value, eaf.outcome = AF_Allele2) %>%
      mutate(id.outcome = disease, outcome = disease) %>%
      select(all_of(c("SNP", "beta.outcome", "se.outcome", "eaf.outcome",
                      "effect_allele.outcome", "other_allele.outcome",
                      "pval.outcome", "id.outcome", "outcome")))
    
    dat <- harmonise_data(exposure_dat = exp, outcome_dat = out, action = 2)
    if (nrow(dat) == 0) { message("    ERROR: No harmonized SNPs."); return(NULL) }
    
    
    mr_res_orig <- mr(dat)
    if (nrow(mr_res_orig) == 0) { message("    ERROR: Original MR analysis failed."); return(NULL) }
    
    mr_res_final <- mr_res_orig
    
    mr_res_final$presso_p <- NA_real_
    
    mr_res_final$nsnp2 <- NA_integer_
    mr_res_final$b2 <- NA_real_
    mr_res_final$se2 <- NA_real_
    mr_res_final$pval2 <- NA_real_
    mr_res_final$presso_p2 <- NA_real_
    
    if (nrow(dat) > 3) { 
      message(paste0("    (Running MR-PRESSO with ", nrow(dat), " SNPs...)"))
      
      presso_res_orig <- NULL
      set.seed(1010)
      presso_res_orig <- tryCatch({
        mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                  SdOutcome = "se.outcome", SdExposure = "se.exposure",
                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                  data = as.data.frame(dat), NbDistribution = 1000,
                  SignifThreshold = 0.05)
      }, error = function(e) {
        message(paste0("    WARNING: MR-PRESSO failed: ", e$message))
        return(NULL)
      })
      
      presso_global_pval_numeric_before <- tryCatch(presso_res_orig$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA_real_)
      
      mr_res_final$presso_p <- presso_global_pval_numeric_before
      
      if (!is.null(presso_res_orig)) {
        presso_path_before <- file.path(final_result_dir, paste0("MR_PRESSO_before_", disease, "_", tissue_file_id, ".rds"))
        saveRDS(presso_res_orig, presso_path_before)
      }
      
      outlier_idx <- tryCatch(presso_res_orig$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, error = function(e) NULL)
      
      if (!is.na(presso_global_pval_numeric_before) && presso_global_pval_numeric_before < 0.05 && 
          !is.null(outlier_idx) && is.numeric(outlier_idx) && length(outlier_idx) > 0) {
        
        message(paste("    Found", length(outlier_idx), "outliers (PRESSO p < 0.05). Running corrected MR analysis..."))
        dat_no_outliers <- dat[-outlier_idx, ]
        
        presso_global_pval_numeric_after <- NA_real_
        
        if (nrow(dat_no_outliers) > 2) { 
          mr_res_corr <- mr(dat_no_outliers)
          
          presso_res_after <- NULL
          if (nrow(dat_no_outliers) > 3) { 
            set.seed(1010)
            presso_res_after <- tryCatch({
              mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure",
                        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                        data = as.data.frame(dat_no_outliers), NbDistribution = 1000,
                        SignifThreshold = 0.05)
            }, error = function(e) NULL)
            
            presso_global_pval_numeric_after <- tryCatch(presso_res_after$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA_real_)
          }
          
          if (!is.null(presso_res_after)) {
            presso_path_after <- file.path(final_result_dir, paste0("MR_PRESSO_after_", disease, "_", tissue_file_id, ".rds"))
            saveRDS(presso_res_after, presso_path_after)
          }
          
          mr_res_corr_renamed <- mr_res_corr %>%
            select(id.exposure, id.outcome, outcome, exposure, method,
                   b2 = b, 
                   se2 = se, 
                   pval2 = pval, 
                   nsnp2 = nsnp)
          
          mr_res_final <- full_join(mr_res_orig, mr_res_corr_renamed, 
                                    by = c("id.exposure", "id.outcome", "outcome", "exposure", "method"))
          
          mr_res_final$presso_p <- presso_global_pval_numeric_before
          
          mr_res_final$presso_p2 <- presso_global_pval_numeric_after
          
        } else {
          message("    Insufficient SNPs after outlier removal (<= 2). Corrected MR analysis cannot be performed.")
        }
        
      } else {
        message(paste("    Skipping corrected MR analysis: (PRESSO pval:", 
                      if(is.na(presso_global_pval_numeric_before)) "NA" else round(presso_global_pval_numeric_before, 3), 
                      ", Outliers:", if(is.null(outlier_idx)) 0 else length(outlier_idx), ")"))
      }
      
    } else {
      message(paste0("    (Skipping MR-PRESSO: only ", nrow(dat), " SNP(s) available.)"))
      message(paste("    Skipping corrected MR analysis: (PRESSO pval: NA , Outliers: 0 )"))
    }
    
    mr_res_final$tissue <- tissue
    mr_res_final$disease <- disease
    
    mr_res_final <- mr_res_final %>%
      mutate(
        presso_p = fmt_pval(presso_p),
        presso_p2 = fmt_pval(presso_p2)
      )
    
    out_path <- file.path(final_result_dir, paste0("MRresult_with_PRESSO_", disease, "_", tissue_file_id, ".csv"))
    fwrite(mr_res_final, out_path)
    
    return(mr_res_final) 
    
  }, error = function(e) {
    message(paste("  ERROR in MR run for index", i, ":", e$message))
    return(NULL) # Return NULL on error
  })
}, future.seed = TRUE)

# --- Combine All Results into One File ---
valid_results <- all_results_list[lengths(all_results_list) > 0]
if (length(valid_results) > 0) {
  final_combined <- bind_rows(valid_results)
  final_output_file <- file.path(mr_output_dir, "MR_all_results_corrected_combined.csv")
  fwrite(final_combined, final_output_file)
  message(paste("\n🎉 --- ANALYSIS COMPLETE! --- 🎉\nFinal combined results (with corrected) saved to:", final_output_file))
} else {
  message("\n❌ --- ANALYSIS FAILED --- ❌\nNo MR results were generated.")
}