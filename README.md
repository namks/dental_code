### Overview

This repository contains the analysis scripts used to perform genome-wide association studies (GWAS), gene-based rare variant tests, cross-trait association analyses, Mendelian randomization (MR), LD score regression (LDSC), and adjusted odds ratio estimation.
Each script corresponds to a specific component of the analysis workflow.

⸻

### File Descriptions

#### 1. `gwas_step1.sh`

Performs SAIGE step 1 (null model fitting) for binary traits.
Includes covariate specification, sparse GRM usage, and null GLMM model generation per phenotype.

⸻

#### 2. `gwas_step2.sh`

Runs SAIGE step 2 (single-variant GWAS) using imputed genotype data.
Produces variant-level association statistics across chromosomes.

⸻

#### 3. `gwas_step2_gene.sh`

Executes SAIGE gene-based tests (ExWAS) for rare variants.
Uses WES-based PLINK files and gene-level group files.

⸻

#### 4. `ldsc.sh`

Conducts LD score regression for:
	•	SNP heritability estimation
	•	Genetic correlation between phenotypes

Includes summary statistic munging and LDSC model execution.

⸻

#### 5. `cpassoc.R`

Implements CPASSOC for cross-trait SNP association testing:
	•	Computes S_Hom and S_Het statistics
	•	Combines evidence from two traits
	•	Identifies shared genetic signals

⸻

#### 6. `mr.R`

Runs TwoSampleMR + MR-PRESSO for causal inference across disease pairs:
	•	Harmonization of exposure/outcome GWAS
	•	IVW, MR-Egger, Weighted median, etc.
	•	Global outlier test (MR-PRESSO)
	•	Filtered and unfiltered MR analyses

⸻

#### 7. `adjusted_OR.R`

Computes adjusted odds ratios for dental–systemic disease pairs:
	•	Logistic regression with Firth correction
	•	Adjustment for demographic, clinical, and lifestyle covariates
	•	Outputs OR, SE, and p-values for each disease pair

⸻

#### 8. `twas.R`

Runs FUSION TWAS pipeline for transcriptomic causal inference across multiple tissues

⸻

#### 9. `twas_post_process.R`

Generate .top files for transcriptome-wide significant associations and runs FUSION post_process pipeline for Joint/conditional tests

⸻

#### 10. `twas_coloc_permutation.R`

Generate .top files for transcriptome-wide significant associations and runs FUSION post_process pipeline for Joint/conditional tests


