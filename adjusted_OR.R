## R code to calculate adjusted odds ratios with Firth correction
## for associations between dental and systemic diseases

library(data.table)
library(logistf)

## -----------------------------------------------------------------------------
## 1. User-configurable paths
## -----------------------------------------------------------------------------

# Root directory of the project (modify for your environment)
PROJECT_ROOT <- "/path/to/project/root"

# Input phenotype / covariate files (example: UK Biobank-based)
DENTAL_PHENO_FILE   <- file.path(PROJECT_ROOT, "phenotypes", "dental_pheno.txt.gz")
SYSTEMIC_PHENO_FILE <- file.path(PROJECT_ROOT, "phenotypes", "systemic_pheno.txt.gz")
COVARIATE_FILE      <- file.path(PROJECT_ROOT, "phenotypes", "systemic_covariates.csv")

## -----------------------------------------------------------------------------
## 2. Load data
## -----------------------------------------------------------------------------

covar         <- fread(COVARIATE_FILE)      # should contain ID + covariates
systemic_pheno <- fread(SYSTEMIC_PHENO_FILE)
dental_pheno   <- fread(DENTAL_PHENO_FILE)

## In this example:
## - dental_pheno has column "IID" (individual ID) + dental phenotypes
## - systemic_pheno has column "IID" + systemic phenotypes
## - covar has column "eid" (or similar) + covariates (age, sex, BMI, etc.)

# Merge data by individual ID
all_dat <- merge(dental_pheno, systemic_pheno, by = "IID")
all_dat <- merge(all_dat, covar, by.x = "IID", by.y = "eid")

## -----------------------------------------------------------------------------
## 3. Phenotype definitions
## -----------------------------------------------------------------------------

# Dental phenotypes (can be extended to a vector of multiple traits)
dental_pheno_list <- "list_of_orofacial_diseases"  # e.g., "dental_caries", "periodontitis"

# Systemic phenotypes
systemic_pheno_list <- "list_of_systemic_diseases" # e.g., "cad", "all_cancer"

## -----------------------------------------------------------------------------
## 4. Covariate definitions
##    (example: UK Biobank field IDs are kept as-is but can be relabeled)
## -----------------------------------------------------------------------------
## For clarity, you may want to rename these in a preprocessing step, e.g.
##   sex, age, bmi, smoking, income, alcohol, LDL, HDL, TG, CRP, etc.
## Here we keep the original variable names but group them into lists.

# Dental-side covariates (example for dental_caries)
dental_covar <- list(
  "list_of_dental_disease1_covariates", "list_of_dental_disease2_covariates", ... 
)

# Systemic-side covariates for each systemic phenotype
systemic_covar <- list(
  "list_of_systemic_disease1_covariates", "list_of_systemic_disease2_covariates", ...
)

## -----------------------------------------------------------------------------
## 5. Containers for results
## -----------------------------------------------------------------------------

n_dental   <- length(dental_pheno_list)
n_systemic <- length(systemic_pheno_list)

adjusted_or_df <- data.frame(
  matrix(NA_real_, nrow = n_dental, ncol = n_systemic),
  row.names = dental_pheno_list
)
colnames(adjusted_or_df) <- systemic_pheno_list

pval_df <- adjusted_or_df
se_df   <- adjusted_or_df

## -----------------------------------------------------------------------------
## 6. Fit Firth-corrected logistic regression models
## -----------------------------------------------------------------------------

for (i in seq_along(dental_pheno_list)) {
  for (j in seq_along(systemic_pheno_list)) {
    
    dental_trait   <- dental_pheno_list[i]
    systemic_trait <- systemic_pheno_list[j]
    
    dental_covar_list   <- dental_covar[[i]]
    systemic_covar_list <- systemic_covar[[j]]
    
    # Union of covariates from exposure and outcome side
    covar_list <- unique(c(dental_covar_list, systemic_covar_list))
    
    # Build formula string: systemic ~ dental + covariates
    rhs <- paste(c(dental_trait, covar_list), collapse = " + ")
    fml <- as.formula(paste(systemic_trait, "~", rhs))
    
    # Firth logistic regression
    fit_firth <- logistf(fml, data = all_dat)
    
    # Coefficient for the dental phenotype (2nd term in the model)
    beta_hat <- coef(fit_firth)[dental_trait]
    se_hat   <- sqrt(diag(vcov(fit_firth)))[dental_trait]
    p_hat    <- fit_firth$prob[which(names(coef(fit_firth)) == dental_trait)]
    
    adjusted_or_df[i, j] <- exp(beta_hat)
    pval_df[i, j]        <- p_hat
    se_df[i, j]          <- se_hat
  }
  message("Finished dental phenotype index: ", i)
}

## -----------------------------------------------------------------------------
## 7. Inspect / save results
## -----------------------------------------------------------------------------

print(adjusted_or_df)
print(pval_df)
print(se_df)

# Optionally save to file
OUT_DIR <- file.path(PROJECT_ROOT, "results", "firth_or")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

fwrite(adjusted_or_df, file.path(OUT_DIR, "adjusted_odds_ratio_firth.txt"),
       sep = "\t", quote = FALSE, row.names = TRUE)
fwrite(pval_df,        file.path(OUT_DIR, "pvalues_firth.txt"),
       sep = "\t", quote = FALSE, row.names = TRUE)
fwrite(se_df,          file.path(OUT_DIR, "se_firth.txt"),
       sep = "\t", quote = FALSE, row.names = TRUE)
