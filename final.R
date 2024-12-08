setwd("~/academics/DSC190/project")

library(data.table)
library(tidyr)
library(dplyr)

# read in gene expression file, filter by hsq.pv < 0.05 to cut down the number of genes
significant_hsq_genes <- c()
files <- list.files(path="fusion_twas-master/WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue", pattern="*.wgt.RDat", full.names=TRUE, recursive=FALSE)
for (i in files) {
  load(i)
  if (hsq.pv < 0.05) {
    significant_hsq_genes <- append(significant_hsq_genes, i)
  }
}

# FIX GWAS COLUMN NAMES

# read in gwas 1
gwas <- fread("fusion_twas-master/34737426-GCST90043925-EFO_1000326.h.tsv.gz") %>%
  mutate(Z = hm_beta / standard_error) %>%
  rename(SNP = hm_rsid, A1 = hm_effect_allele, A2 = hm_other_allele) %>%
  mutate(log_10_value = -log(x = p_value, base = 10))

write.table(gwas, "fusion_twas-master/GCST90043925_harmonized.tsv.gz")

plot(gwas$log_10_value)

# read in gwas 2
gwas <- fread("fusion_twas-master/GCST90454347.h.tsv.gz") %>%
  mutate(Z = beta / standard_error) %>%
  rename(SNP = rsid, A1 = effect_allele, A2 = other_allele)

write.table(gwas, "fusion_twas-master/harmonized2.tsv", row.names=F, col.names=T)

plot(gwas$neg_log_10_p_value)
abline(h = -log(x=0.05, base=10), col="red")

# RUN FUSION

# read in TWAS results: gwas 1 & GTEx eQTL
twasfiles <- list.files(path="fusion_twas-master", pattern="*_gtex.dat", full.names=TRUE, recursive=FALSE)
twas_significant_genes <- matrix(0, 0, ncol(twaschr1) + 1)
colnames(twas_significant_genes) <- c(colnames(twaschr1),"FDR")
for (i in twasfiles) {
  twaschri <- fread(i)
  twaschri$FDR <- p.adjust(twaschri$TWAS.P, method="BH")
  twas_significant_genes <- rbind(twas_significant_genes, twaschri[which(twaschri$FDR < 0.05),])
  plot(twaschri$FDR, main=i)
  abline(h = 0.05, col="red")
}

write.table(twas_significant_genes, "fusion_twas-master/breastcancer_gtex.top")

# read in TWAS results: gwas 1 & TCGA eQTL
twasfiles <- list.files(path="fusion_twas-master", pattern="*_tcga_tumor_exp.dat", full.names=TRUE, recursive=FALSE)
twaschr1 <- fread(twasfiles[1])
twas_significant_genes <- matrix(0, 0, ncol(twaschr1) + 1)
colnames(twas_significant_genes) <- c(colnames(twaschr1),"FDR")
for (i in twasfiles) {
  twaschri <- fread(i)
  twaschri$FDR <- p.adjust(twaschri$TWAS.P, method="BH")
  twas_significant_genes <- rbind(twas_significant_genes, twaschri[which(twaschri$FDR < 0.05),])
  plot(twaschri$FDR, main=i)
  abline(h = 0.05, col="red")
}

write.table(twas_significant_genes, "fusion_twas-master/breastcancer_tcga_tumor_ex.top")

# read in TWAS results: gwas 2 & GTEx eQTL
twasfiles <- list.files(path="fusion_twas-master", pattern="*_gtex2.dat", full.names=TRUE, recursive=FALSE)
twaschr1 <- fread(twasfiles[1])
twas_significant_genes <- matrix(0, 0, ncol(twaschr1) + 1)
colnames(twas_significant_genes) <- c(colnames(twaschr1),"FDR")
for (i in twasfiles) {
  twaschri <- fread(i)
  twaschri$FDR <- p.adjust(twaschri$TWAS.P, method="BH")
  twas_significant_genes <- rbind(twas_significant_genes, twaschri[which(twaschri$FDR < 0.05),])
  plot(twaschri$FDR, main=i)
  abline(h = 0.05, col="red")
}

# read in TWAS results: gwas 2 & TCGA eQTL
twasfiles <- list.files(path="fusion_twas-master", pattern="*_tcga_tumor_exp2.dat", full.names=TRUE, recursive=FALSE)
twaschr1 <- fread(twasfiles[1])
twas_significant_genes <- matrix(0, 0, ncol(twaschr1) + 1)
colnames(twas_significant_genes) <- c(colnames(twaschr1),"FDR")
for (i in twasfiles) {
  twaschri <- fread(i)
  twaschri$FDR <- p.adjust(twaschri$TWAS.P, method="BH")
  twas_significant_genes <- rbind(twas_significant_genes, twaschri[which(twaschri$FDR < 0.05),])
  plot(twaschri$FDR, main=i)
  abline(h = 0.05, col="red")
}

# FOCUS fine mapping
# match SNPs for each chromosome
twas_sig_locus_chr <- c(1, 2, 5, 7, 11, 12, 15, 17, 19, 22)
gwas1 <- fread("fusion_twas-master/GCST90043925_harmonized.tsv.gz", header=F)
gwas2 <- fread("fusion_twas-master/harmonized2.tsv", header=F)
for (i in twas_sig_locus_chr) {
  gwas1_i <- copy(gwas1)
  gwas2_i <- copy(gwas2)
  snps <- fread(sprintf("fusion_twas-master/LDREF/1000G.EUR.%s.bim", i), header=F)
  w1 <- which(gwas1_i$SNP %in% snps$V2)
  w2 <- which(gwas2_i$SNP %in% snps$V2)
  
  gwas1_i <- gwas1_i[w1, ]
  gwas2_i <- gwas2_i[w2, ]
  
  gwas1_i$CHR <- i
  gwas2_i$CHR <- i
  gwas1_i$BP <- snps$V4[match(gwas1_i$SNP, snps$V2)]
  gwas2_i$BP <- snps$V4[match(gwas2_i$SNP, snps$V2)]
  
  write.table(gwas1_i, file=sprintf("harmonized1.sumstats.chr%s", i), row.names=F, col.names=T, sep="\t", quote=F)
  write.table(gwas2_i, file=sprintf("harmonized2.sumstats.chr%s", i), row.names=F, col.names=T, sep="\t", quote=F)
}

# COLOC
library(coloc)
gtex_significant_genes <- fread("fusion_twas-master/breastcancer_gtex.top")

# Gene: ENSG00000255274.9 (Chr 11)
w_SMIM35 <- list.files(path="fusion_twas-master/WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue", pattern="ENSG00000255274.9*", full.names=TRUE, recursive=FALSE)

D1_eqtl <- list()

# FOCUS fine mapping
# chr 2
y <- fread("fusion_twas-master/GCST90043925_harmonized.tsv.gz", header=T)
y <- fread("fusion_twas-master/GCST90454347.v2.tsv.gz", header=T)
y <- fread("fusion_twas-master/harmonized2.tsv", header=T)
y <- y[, 2:20]
colnames(y) <- colnames(gwas)

write.table(y, "fusion_twas-master/harmonized2.tsv")
snps <- fread("fusion_twas-master/LDREF/1000G.EUR.2.bim", header=F)
# w <- which(y$hm_rsid %in% snps$V2)
w <- which(y$SNP %in% snps$V2)
y <- y[w,]
y$CHR <- 2
y$BP <- snps$V4[match(y$SNP, snps$V2)]
y$SNP <- y$SNP
y$A1 <- y$hm_effect_allele
y$A2 <- y$hm_other_allele
y$Z <- y$hm_beta / y$standard_error
write.table(y, file="harmonized2.sumstats.chr2", row.names=F, col.names=T, sep="\t", quote=F)

# chr 7
y <- fread("fusion_twas-master/34737426-GCST90043925-EFO_1000326.h.tsv.gz", header=T)
snps <- fread("fusion_twas-master/LDREF/1000G.EUR.7.bim", header=F)
w <- which(y$hm_rsid %in% snps$V2)
y <- y[w,]
y$CHR <- 7
y$BP <- snps$V4[match(y$hm_rsid, snps$V2)]
y$SNP <- y$hm_rsid
y$A1 <- y$hm_effect_allele
y$A2 <- y$hm_other_allele
y$Z <- y$hm_beta / y$standard_error
write.table(y, file="harmonized.sumstats.chr7", row.names=F, col.names=T, sep="\t", quote=F)

# chr 11
y <- fread("fusion_twas-master/34737426-GCST90043925-EFO_1000326.h.tsv.gz", header=T)
snps <- fread("fusion_twas-master/LDREF/1000G.EUR.11.bim", header=F)
w <- which(y$hm_rsid %in% snps$V2)
y <- y[w,]
y$CHR <- 11
y$BP <- snps$V4[match(y$hm_rsid, snps$V2)]
y$SNP <- y$hm_rsid
y$A1 <- y$hm_effect_allele
y$A2 <- y$hm_other_allele
y$Z <- y$hm_beta / y$standard_error
write.table(y, file="harmonized.sumstats.chr11", row.names=F, col.names=T, sep="\t", quote=F)

load("fusion_breastmammary.db")

# TCGA Tumor Expression Predictive Model


# chr 5
y <- fread("fusion_twas-master/34737426-GCST90043925-EFO_1000326.h.tsv.gz", header=T)
snps <- fread("fusion_twas-master/LDREF/1000G.EUR.5.bim", header=F)
w <- which(y$hm_rsid %in% snps$V2)
y <- y[w,]
y$CHR <- 5
y$BP <- snps$V4[match(y$hm_rsid, snps$V2)]
y$SNP <- y$hm_rsid
y$A1 <- y$hm_effect_allele
y$A2 <- y$hm_other_allele
y$Z <- y$hm_beta / y$standard_error
write.table(y, file="harmonized.sumstats.chr5", row.names=F, col.names=T, sep="\t", quote=F)

# chr 22
y <- fread("fusion_twas-master/34737426-GCST90043925-EFO_1000326.h.tsv.gz", header=T)
snps <- fread("fusion_twas-master/LDREF/1000G.EUR.22.bim", header=F)
w <- which(y$hm_rsid %in% snps$V2)
y <- y[w,]
y$CHR <- 22
y$BP <- snps$V4[match(y$hm_rsid, snps$V2)]
y$SNP <- y$hm_rsid
y$A1 <- y$hm_effect_allele
y$A2 <- y$hm_other_allele
y$Z <- y$hm_beta / y$standard_error
write.table(y, file="harmonized.sumstats.chr22", row.names=F, col.names=T, sep="\t", quote=F)
