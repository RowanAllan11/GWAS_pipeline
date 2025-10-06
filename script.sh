#!/bin/bash

# GWAS Pipeline
# Author: Rowan Allan
# Date: 01/11/24
# Dependencies: plink (v1.9+), R (with qqman)
# Input folder: input/ (gwas.bed, gwas.bim, gwas.fam, gwas.covar)

set -e  # exit on error

# 1. Initial Data Overview
plink --bfile input/gwas --freq --out summary
awk '{print $6}' input/gwas.fam | sort | uniq -c  # count cases/controls

# 2. Missing Data QC
plink --bfile input/gwas --missing --out missing.stat
grep A2038 miss.stat.imiss          # missing SNPs for individual A2038
grep rs2493272 miss.stat.lmiss      # missing individuals for SNP rs2493272

# 3. Minor Allele Frequencies
grep rs4970357 summary.frq

# 4. Apply QC Filters
plink --bfile input/gwas --geno 0.05 --make-bed --out SNP.miss
plink --bfile SNP.miss --maf 0.01 --make-bed --out QC.data
plink --bfile QC.data --missing --out QC.vis

# 5. Basic Association Testing
plink --bfile QC.data --snp rs9651273 --model --out model.results

# 6. Association Testing with Multiple Correction
plink --bfile QC.data --assoc --adjust --out assoc.adjusted

# 7. Logistic Regression with Covariates
cp input/gwas.covar gwas.covar.txt
plink --bfile QC.data --logistic --covar gwas.covar.txt --sex --adjust --out yes

# 8. R Analysis & Visualization
Rscript analysis.R
