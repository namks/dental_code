## Example R script to run two-sample MR and MR-PRESSO
## using pre-formatted exposure and outcome summary statistics

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MRPRESSO)

options(width = 200)

## -----------------------------
## Phenotype definitions
## -----------------------------

# Exposure phenotypes
exp_pheno <- "list_of_exposure_phenotypes"

# Outcome phenotypes
out_pheno <- "list_of_outcome_phenotypes"

## -----------------------------
## Main loop over exposure/outcome pairs
## -----------------------------

for (ex in exp_pheno) {
  for (ou in out_pheno) {
    
    message(sprintf("Running MR: %s -> %s", ex, ou))
    
    try({
      ## -----------------------------
      ## Load exposure & outcome data
      ## -----------------------------
      
      f_exp <- file.path(MR_INPUT_DIR, paste0(ex, ".txt"))
      f_out <- file.path(MR_INPUT_DIR, paste0(ou, ".txt"))
      
      exposure <- fread(f_exp)  |> as.data.frame()
      outcome  <- fread(f_out)  |> as.data.frame()
      
      ## -----------------------------
      ## Format exposure data
      ## -----------------------------
      
      exposure$pheno <- ex
      exposure$N     <- exposure$N_case + exposure$N_ctrl
      
      exposure_dat <- format_data(
        type                = "exposure",
        phenotype_col       = "pheno",
        dat                 = exposure,
        snps                = NULL,
        snp_col             = "SNP",
        beta_col            = "BETA",
        se_col              = "SE",
        pval_col            = "p.value",
        effect_allele_col   = "Allele2",
        other_allele_col    = "Allele1",
        eaf_col             = "AF_Allele2",
        samplesize_col      = "N"
      )
      
      ## -----------------------------
      ## Format outcome data
      ## -----------------------------
      
      outcome$pheno <- ou
      outcome$N     <- outcome$N_case + outcome$N_ctrl
      
      outcome_dat <- format_data(
        type                = "outcome",
        phenotype_col       = "pheno",
        dat                 = outcome,
        # snps              = exposure$SNP,  # optional restriction
        snp_col             = "SNP",
        beta_col            = "BETA",
        se_col              = "SE",
        pval_col            = "p.value",
        effect_allele_col   = "Allele2",
        other_allele_col    = "Allele1",
        eaf_col             = "AF_Allele2",
        samplesize_col      = "N"
      )
      
      ## -----------------------------
      ## Harmonisation
      ## -----------------------------
      
      dat <- harmonise_data(
        exposure_dat = exposure_dat,
        outcome_dat  = outcome_dat
      )
      
      ## -----------------------------
      ## MR (no filtering)
      ## -----------------------------
      
      res_no_filter <- mr(dat)
      
      ## MR-PRESSO (no filtering)
      t1_no_filter <- mr_presso(
        BetaOutcome    = "beta.outcome",
        BetaExposure   = "beta.exposure",
        SdOutcome      = "se.outcome",
        SdExposure     = "se.exposure",
        OUTLIERtest    = TRUE,
        DISTORTIONtest = TRUE,
        data           = dat,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
      
      mr_presso_pval_no_filter <-
        t1_no_filter$`MR-PRESSO results`$`Global Test`$Pvalue
      
      print(t1_no_filter)
      
      ## -----------------------------
      ## MR after filtering on outcome p-value
      ## (e.g. remove IVs strongly associated with outcome)
      ## -----------------------------
      
      dat_filt <- dat[dat$pval.outcome > (0.05 / nrow(dat)), ]
      
      res_filter <- mr(dat_filt)
      
      t1_filter <- mr_presso(
        BetaOutcome    = "beta.outcome",
        BetaExposure   = "beta.exposure",
        SdOutcome      = "se.outcome",
        SdExposure     = "se.exposure",
        OUTLIERtest    = TRUE,
        DISTORTIONtest = TRUE,
        data           = dat_filt,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
      
      mr_presso_pval_filter <-
        t1_filter$`MR-PRESSO results`$`Global Test`$Pvalue
      
      print(t1_filter)
      
      ## -----------------------------
      ## Save results
      ## -----------------------------
      
      outname_no_filter <- file.path(
        MR_OUTPUT_DIR,
        paste0(ex, "_", ou, ".txt")
      )
      outname_filter <- file.path(
        MR_OUTPUT_DIR,
        paste0(ex, "_", ou, "_filtered.txt")
      )
      outname_mr_presso_no_filter <- file.path(
        MR_OUTPUT_DIR,
        paste0(ex, "_", ou, "_mr_presso.txt")
      )
      outname_mr_presso_filter <- file.path(
        MR_OUTPUT_DIR,
        paste0(ex, "_", ou, "_filtered_mr_presso.txt")
      )
      
      write.table(res_no_filter, outname_no_filter,
                  sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(res_filter, outname_filter,
                  sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(mr_presso_pval_no_filter, outname_mr_presso_no_filter,
                  sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(mr_presso_pval_filter, outname_mr_presso_filter,
                  sep = "\t", quote = FALSE, row.names = FALSE)
    })
  }
}
