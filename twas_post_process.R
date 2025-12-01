################################################################################
# R code to perform TWAS post process
#
# Steps:
# 3. Generate .top files
# 4. Post-Process
# 5. Merge results
################################################################################

suppressMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
})

cfg <- list(
  twas_output_dir = "./results/",
  sumstats_split_dir = "./SUMSTATS/split/",
  post_input_top_dir = "./results/top_files/",
  fusion_post_script = "./fusion_twas-master/FUSION.post_process.R",
  ld_ref_prefix = "./fusion_twas-master/LDREF/1000G.EUR.",
  post_merged_output_file = "./results/top_files/all_joint_included_results_for_coloc.txt"
)

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
  flush.console()
}

# ---------------------------------------------------------------------------- #
# Step 3: Generate .top Files
# ---------------------------------------------------------------------------- #
generate_top_files <- function(cfg) {
  log_message("--- [Generating .top files from TWAS Results] ---")
  
  if (!dir.exists(cfg$post_input_top_dir)) dir.create(cfg$post_input_top_dir, recursive = TRUE)
  
  filtered_files <- list.files(cfg$twas_output_dir, pattern = "_All_Tissues_TWAS_Filtered_Results.txt$", full.names = TRUE)
  
  if(length(filtered_files) == 0) stop("No filtered results found from TWAS.")
  
  count <- 0
  for (f in filtered_files) {
    dt <- fread(f)
    disease_name <- gsub("_All_Tissues_TWAS_Filtered_Results.txt", "", basename(f))
    setnames(dt, names(dt), gsub(" ", "_", names(dt)))
    
    char_cols <- names(dt)[sapply(dt, is.character)]
    for (col in char_cols) {

      dt[[col]] <- gsub(" ", "_", dt[[col]])

      dt[get(col) == "", (col) := "NA"]

      dt[is.na(get(col)), (col) := "NA"]
    }
    
    num_cols <- names(dt)[sapply(dt, is.numeric)]
    for (col in num_cols) {
    }
    # ---------------------------
    
    groups <- split(dt, by = c("Tissue", "CHR"))
    
    for (grp_name in names(groups)) {
      sub_dt <- groups[[grp_name]]
      if (nrow(sub_dt) == 0) next
      
      tissue <- sub_dt$Tissue[1]
      chr <- sub_dt$CHR[1]
      
      top_filename <- paste0(disease_name, ".", tissue, ".chr", chr, ".dat.top")
      top_path <- file.path(cfg$post_input_top_dir, top_filename)
      
      fwrite(sub_dt, top_path, sep = "\t", quote = FALSE, na = "NA")
      count <- count + 1
    }
  }
  log_message(sprintf("Generated %d .top files.", count))
}

# ---------------------------------------------------------------------------- #
# Step 4: Post-Process (Conditional Analysis)
# ---------------------------------------------------------------------------- #
run_step_4_post_process <- function(cfg) {
  log_message("--- [Step 4: Running FUSION Post-Process] ---")
  
  top_files <- list.files(cfg$post_input_top_dir, pattern = "\\.top$", full.names = TRUE)
  
  if (length(top_files) == 0) {
    log_message("No .top files found. Check generation step.")
    return()
  }
  
  for (top_file in top_files) {
    top_base <- basename(top_file)
    
    parts <- strsplit(top_base, "\\.")[[1]]
    disease_name <- parts[1]
    
    chr_idx <- grep("^chr\\d+", parts)[1]
    tissue_name <- paste(parts[2:(chr_idx-1)], collapse = ".")
    chr <- as.integer(sub("chr", "", parts[chr_idx]))
    
    sumstats_file <- file.path(cfg$sumstats_split_dir, paste0(disease_name, "_chr", chr, ".sumstats"))
    
    if (!file.exists(sumstats_file)) next
    
    out_prefix <- sub("\\.top$", ".top.analysis", top_file)
    included_file <- paste0(out_prefix, ".joint_included.dat")
    
    if (file.exists(included_file)) next
    
    cmd <- paste(
      "Rscript", shQuote(cfg$fusion_post_script),
      "--sumstats", shQuote(sumstats_file),
      "--input", shQuote(top_file),
      "--out", shQuote(out_prefix),
      "--ref_ld_chr", shQuote(cfg$ld_ref_prefix),
      "--chr", chr,
      "--plot --locus_win 100000"
    )
    
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    cat(".") 
  }
  cat("\n")
  log_message("--- [Step 4 Complete] ---")
}

# ---------------------------------------------------------------------------- #
# Step 5: Merge Results & Prepare for Step 6
# ---------------------------------------------------------------------------- #
run_step_5_merge_post_process <- function(cfg) {
  log_message("--- [Step 5: Merging Joint Included Results] ---")
  
  files <- list.files(cfg$post_input_top_dir, pattern = "\\.joint_included\\.dat$", full.names = TRUE)
  
  if (length(files) == 0) {
    log_message("No joint_included files found.")
    return()
  }
  
  joint_list <- list()
  
  for (f in files) {
    dt <- fread(f)
    fname <- basename(f)
    
    parts <- strsplit(fname, "\\.")[[1]]
    
    chr_idx <- grep("^chr\\d+", parts)[1]
    disease_name <- parts[1]
    tissue_name <- paste(parts[2:(chr_idx-1)], collapse = ".")
    chr_val <- as.integer(sub("chr", "", parts[chr_idx]))
    
    dt[, Disease := disease_name]
    dt[, Tissue := tissue_name]
    dt[, CHR := chr_val]
    
    joint_list[[length(joint_list) + 1]] <- dt
  }
  
  joint_all <- rbindlist(joint_list, fill = TRUE)
  
  fwrite(joint_all, cfg$post_merged_output_file, sep = "\t")
  log_message(paste("✅ Saved Joint Significant Results to:", cfg$post_merged_output_file))
}

# Execution
if (!interactive()) {
  generate_top_files(cfg)
  run_step_4_post_process(cfg)
  run_step_5_merge_post_process(cfg)

}


