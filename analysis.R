# analysis.R
# R Script for GWAS QC, Visualization, and Significant SNP Identification
# Author: Rowan Allan
# Date: 01/11/24

# Load required packages
library(qqman)

# -------------------------------
# 1. Missingness Plots
# -------------------------------
lmiss_data <- read.table("miss.stat.lmiss", header = TRUE)
imiss_data <- read.table("miss.stat.imiss", header = TRUE)

# SNP missingness histogram
png("plots/snp_missingness.png", width=800, height=600)
hist(lmiss_data$F_MISS,
     main = "Histogram of SNP Missingness",
     xlab = "Fraction of Missing Genotypes",
     col = "skyblue",
     breaks = 50)
dev.off()

# Individual missingness histogram
png("plots/ind_missingness.png", width=800, height=600)
hist(imiss_data$F_MISS,
     main = "Histogram of Individual Missingness",
     xlab = "Fraction of Missing Genotypes",
     col = "lightcoral",
     breaks = 30)
dev.off()

# -------------------------------
# 2. Minor Allele Frequency (MAF) Distribution
# -------------------------------
allele_freq <- read.table("summary.frq", header = TRUE)

png("plots/maf_distribution.png", width=800, height=600)
hist(allele_freq$MAF,
     main = "Distribution of Minor Allele Frequencies (MAF)",
     xlab = "Minor Allele Frequency",
     col = "skyblue",
     breaks = 50,
     border = "white")
dev.off()

# -------------------------------
# 3. Identify Significant SNPs
# -------------------------------
assoc_data <- read.table("assoc.adjusted.assoc.adjusted", header = TRUE)

# Example filtering for multiple correction significance (BONF, HOLM, SIDAK, FDR)
significant_snps <- subset(assoc_data,
                           BONF < 0.05 &
                           HOLM < 0.05 &
                           SIDAK_SS < 0.05 &
                           SIDAK_SD < 0.05 &
                           FDR_BH < 0.05 &
                           FDR_BY < 0.05)

write.table(significant_snps, "plots/significant_snps.txt", row.names = FALSE, quote = FALSE)

# -------------------------------
# 4. Manhattan Plot from logistic regression
# -------------------------------
logistic_results <- read.table("yes.assoc.logistic", header = TRUE)

png("plots/manhattan_plot.png", width=1200, height=600)
manhattan(logistic_results,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          main = "Manhattan Plot of Logistic Regression Results")
dev.off()

# -------------------------------
# Finished
# -------------------------------
cat("R analysis complete. Plots saved in 'plots/' folder and significant SNPs listed in 'significant_snps.txt'.\n")
