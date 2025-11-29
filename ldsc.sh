python ${LDSC_DIR}/munge_sumstats.py \
  --sumstats ${PROJECT_ROOT}/results/step2/step2_${pheno[${j}]}_chr${i}_EUR.txt \
  --snp MarkerID \
  --a1 Allele1 \
  --a2 Allele2 \
  --N ${sample[${j}]} \
  --p p.value \
  --out ${PROJECT_ROOT}/results/LDSC/summary/${pheno[${j}]}_chr${i} \
  --chunksize 5000 \
  --merge-alleles ${HM3_SNPLIST}

python ${LDSC_DIR}/ldsc.py \
  --rg \
    ${PROJECT_ROOT}/results/LDSC/summary/${d}_all.sumstats.gz,\
    ${PROJECT_ROOT}/results/LDSC/summary/${s}_all.sumstats.gz \
  --ref-ld-chr ${REF_LD_DIR}/ \
  --w-ld-chr  ${W_LD_DIR}/ \
  --out ${PROJECT_ROOT}/results/LDSC/gc/${d}_${s}
