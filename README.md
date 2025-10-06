# GWAS Pipeline for Association Testing and QC Analysis

## Overview

This project implements a reproducible pipeline for performing Quality Control (QC) and association testing on a GWAS dataset. The analysis utilises **PLINK** for efficient data manipulation and association tests and **R** (specifically the `qqman` package) for visualisation and final result filtering.

The dataset used:

  * **Variants:** 306,102
  * **Individuals:** 4000 (2000 males, 2000 females)
  * **Phenotype Balance:** 2000 Cases and 2000 Controls (balanced cohort)
  * **Population Structure:** 0 non-founders (unrelated cohort)

## Dependencies

The pipeline requires specific tools and R packages defined in `environment.yml`.

| Tool/Package | Version | Purpose |
| :--- | :--- | :--- |
| **plink** | v1.9 | Core tool for genetic data manipulation and association testing. |
| **r-base** | v4.3 | Base R environment. |
| **r-qqman** | Required | R package for generating **Manhattan** and **Q-Q plots**. |


## Setup

```bash
conda env create -f environment.yml
conda activate gwas_pipeline
```

## Usage

The entire workflow is automated via a shell script:

```bash
./run_gwas.sh
```

## Input Data

The pipeline operates on standard PLINK binary files and a separate covariate file. Files not included*

| Filename | Description | Details |
| :--- | :--- | :--- |
| `gwas.bed` | **Binary Genotype Data.** | Requires PLINK to handle. |
| `gwas.bim` | **Variant/SNP Information.** | Includes Chromosome, SNP ID (e.g., rs3934834), base pair position, and alleles (e.g., T/C). |
| `gwas.fam` | **Sample Information.** | Includes Family ID (0), Individual ID (e.g., A2001), Sex (1=male, 2=female), and **Phenotype** (2=case, 1=control). |
| `gwas.covar` | **Covariate Data.** | Contains the **age of the individuals** for adjustment in logistic regression. |


### A. Initial Overview and Quality Control (QC) Checks

| Step | PLINK Command Purpose | Summary of Findings |
| :--- | :--- | :--- |
| **Overview & Frequencies** | Calculate allele frequencies (`--freq`) and confirm case/control balance. | Cohort is balanced (2000 cases/controls). Allele frequencies are calculated. |
| **Missing Data Check** | Generate reports on missingness per SNP (`.lmiss`) and per individual (`.imiss`). | Initial genotyping rate is high at **98.3%**. Missing data examples: Individual A2038 (1.7% missingness) and SNP rs2493272 (2.8% missingness). |
| **Allele Frequencies** | Verify Minor Allele Frequencies (MAF). | Example: SNP rs4970357 has a MAF of 0.05028. Most SNPs have relatively common alleles. |

### B. Filtering and Data Refinement

**Apply QC Filters:** A stringent filtering process is applied to create the cleaned dataset, `QC.data`.

| Filter Type | PLINK Parameter | Threshold | Variants Removed | Rationale |
| :--- | :--- | :--- | :--- | :--- |
| **SNP Missingness** | `--geno` | 0.05 (5%) | **5,552** | Removes poorly genotyped markers. |
| **Minor Allele Frequency** | `--maf` | 0.01 (1%) | **0** | Removes rare variants that provide low statistical power. The high quality of the initial data meant no variants were removed here. |

*Result:* These steps significantly reduced overall missingness and prepared a high-quality variant set for association testing.

### C. Association Testing

**Basic Association Testing (`--model`):**

  * **Test:** Five genetic models were tested for SNP rs9651273 (Additive, Dominant, Recessive, Genotypic, Allelic).
  * *Result:* rs9651273 showed the smallest p-value (**$P=0.01407$**) under the **Dominant (DOM) model**.

**Unadjusted Association Testing and Multiple Correction:**

  * **Test:** General association test (`--assoc`) followed by multiple testing corrections (`--adjust`).
  * *Result:* **9 SNPs** were identified as significant ($P < 0.05$) under all applied corrections (Bonferroni, HOLM, SIDAK, and FDR).
  * **Population Structure Check (Unadjusted):** The Genomic Inflation Factor ($\lambda$) was estimated at **1.00943**. This value is very close to 1, indicating **minimal population stratification** confounding the basic analysis.

**Logistic Regression with Covariates:**

  * **Test:** Logistic regression is performed, adjusting for **sex** and **age** (`--covar` is used after converting `gwas.covar` to `gwas.covar.txt`).
  * *Result:* The $\lambda$ was **$1.0108$**. This value is nearly identical to the unadjusted $\lambda$ ($1.01077$), confirming that **age and sex do not significantly influence the observed population structure or drive the primary genetic associations**.

### D. R Analysis and Visualization

The `analysis.R` script processes the final results from Step 7 to generate plots and define the final list of associated SNPs.

## Expected Outputs and Results

All visualisation and summary files are saved into the `plots/` directory.

| Output File | Content and Significance |
| :--- | :--- |
| `plots/snp_missingness.png` | Histogram showing SNP missingness before QC. |
| `plots/ind_missingness.png` | Histogram showing individual missingness. |
| `plots/maf_distribution.png` | Histogram reflecting the distribution of Minor Allele Frequencies (post-QC). |
| `plots/manhattan_plot.png` | **Manhattan Plot** visualizing the results of the final logistic regression (Step 7). Significant peaks are highlighted on **Chromosomes 3 and 10**.  |
| `plots/qq_plot.png` | **Q-Q Plot** confirming the genomic inflation factor and $P$-value distribution. |
| `plots/significant_snps.txt` | Final list of SNPs that passed the multiple testing thresholds. Examples include **rs6802898, rs7901695, rs7903146, rs8050136**, which are known to be associated with **Type 2 Diabetes (T2D)**. |