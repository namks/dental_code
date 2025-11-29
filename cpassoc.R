## R code to perform CPASSOC to assess cross-trait association

library(data.table)

## -----------------------------
## User-configurable paths
## -----------------------------
# Base project directory (change to your own)
PROJECT_ROOT <- "/path/to/project"

# Directory containing GWAS step2 results
DENTAL_STEP2_DIR   <- file.path(PROJECT_ROOT, "results", "step2", "dental")
SYSTEMIC_STEP2_DIR <- file.path(PROJECT_ROOT, "results", "step2", "systemic")

# Directory to save CPASSOC results
CPASSOC_OUT_DIR <- file.path(PROJECT_ROOT, "results", "CPASSOC")

# Directory where CPASSOC R functions (FunctionSet.R) are stored
CPASSOC_DIR <- "/path/to/CPASSOC"  # contains FunctionSet.R

# Load CPASSOC functions
source(file.path(CPASSOC_DIR, "FunctionSet.R"))

## -----------------------------
## Analysis settings
## -----------------------------

# Exposure phenotype (e.g. dental caries)
exposure <- "exposure_phenotype_name"

# Outcome phenotypes
pheno <- "list_of_systemic_diseases"

for (p in pheno) {
  d <- NULL   # merged container
  
  for (i in 1:22) {
    # Step 2 GWAS results for exposure (e.g. dental caries)
    f_exp <- file.path(
      DENTAL_STEP2_DIR,
      sprintf("step2_%s_chr%d_EUR.txt", exposure, i)
    )
    
    # Step 2 GWAS results for outcome (systemic disease)
    f_out <- file.path(
      SYSTEMIC_STEP2_DIR,
      sprintf("step2_%s_chr%d_EUR.txt", p, i)
    )
    
    d1 <- fread(f_exp)  # exposure
    d2 <- fread(f_out)  # outcome
    
    # Filter by p-value
    d1 <- d1[p.value < 1e-3, ]
    d2 <- d2[p.value < 1e-3, ]
    
    # Merge on SNP ID
    d3 <- merge(d1, d2, by = "MarkerID")
    
    d  <- rbind(d, d3)
    message(sprintf("Processed %s, chr %d", p, i))
  }
  
  ## -----------------------------
  ## Compute Z-scores and CPASSOC
  ## -----------------------------
  
  # Z-scores from two-sided p-values
  d$Z.x <- qnorm(d$p.value.x / 2, lower.tail = FALSE)
  d$Z.y <- qnorm(d$p.value.y / 2, lower.tail = FALSE)
  
  # Assign sign according to beta
  d$Z.x <- ifelse(d$BETA.x > 0,  d$Z.x, -d$Z.x)
  d$Z.y <- ifelse(d$BETA.y > 0,  d$Z.y, -d$Z.y)
  
  # Matrix of Z-scores
  X <- d[, .(Z.x, Z.y)]
  X <- X[complete.cases(X), ]
  
  # Sample sizes (assumes constant per trait)
  N.x <- d$N_case.x[1] + d$N_ctrl.x[1]
  N.y <- d$N_case.y[1] + d$N_ctrl.y[1]
  N   <- c(N.x, N.y)
  
  # Correlation matrix between traits
  R <- cor(X)
  
  # Run CPASSOC SHom
  SHom   <- SHom(X, N, R)
  p_Hom  <- pchisq(SHom, df = 1, ncp = 0, lower.tail = FALSE)
  n_sig_Hom <- sum(p_Hom < 5e-8, na.rm = TRUE)
  message(sprintf("[%s] SHom genome-wide significant loci: %d", p, n_sig_Hom))
  
  # Estimate gamma parameters for SHet
  para <- EstimateGamma(
    N            = 1e4,
    SampleSize   = N,
    CorrMatrix   = R,
    correct      = 1,
    startCutoff  = 0,
    endCutoff    = 1,
    CutoffStep   = 0.05,
    isAllpossible = TRUE
  )
  
  SHet  <- SHet(
    X            = X,
    SampleSize   = N,
    CorrMatrix   = R,
    correct      = 1,
    isAllpossible = TRUE
  )
  
  p_Het <- pgamma(
    q     = SHet - para[3],
    shape = para[1],
    scale = para[2],
    lower.tail = FALSE
  )
  n_sig_Het <- sum(p_Het < 5e-8, na.rm = TRUE)
  message(sprintf("[%s] SHet genome-wide significant loci: %d", p, n_sig_Het))
  
  ## -----------------------------
  ## Save results
  ## -----------------------------
  
  out <- cbind(d, p_Hom = p_Hom, p_Het = p_Het)
  
  if (!dir.exists(CPASSOC_OUT_DIR)) {
    dir.create(CPASSOC_OUT_DIR, recursive = TRUE)
  }
  
  outname <- file.path(
    CPASSOC_OUT_DIR,
    sprintf("CPASSOC_%s_%s.txt", exposure, p)
  )
  
  fwrite(out, outname, sep = "\t", quote = FALSE, row.names = FALSE)
  system(paste("gzip -f", shQuote(outname)))
}
