docker run \
  -v ${PROJECT_ROOT}:${PROJECT_ROOT} \
  wzhou88/saige:1.3.0 step1_fitNULLGLMM.R \
  --sparseGRMFile=${PROJECT_ROOT}/GRM/UKB_sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=${PROJECT_ROOT}/GRM/UKB_sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --plinkFile=${PROJECT_ROOT}/genotype/pruned/prune_all \
  --phenoFile=${PHENO_FILE} \
  --skipVarianceRatioEstimation=FALSE \
  --phenoCol=${TRAIT} \
  --covarColList=Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --LOCO=FALSE \
  --nThreads=16 \
  --outputPrefix=${PROJECT_ROOT}/results/step1/step1_${TRAIT}_EUR \
  --IsOverwriteVarianceRatioFile=TRUE
