library(data.table)
library(tidyr)
library(dplyr)

##### FIX GWAS COLUMN NAMES

gwas <- fread("fusion_twas-master/GCST90454347.h.tsv.gz") %>%
  mutate(Z = beta / standard_error) %>%
  rename(SNP = rsid, A1 = effect_allele, A2 = other_allele)

write.table(gwas, "fusion_twas-master/harmonized2.tsv", row.names=F, col.names=T, sep="\t", quote=F)

##### RUN FUSION (see final.sh)

##### AGGREGATE TWAS FUSION RESULTS ACROSS CHROMOSOMES

# read in TWAS results: GTEx eQTL
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

twas_significant_genes <- twas_significant_genes[order(CHR)]

write.table(twas_significant_genes, "fusion_twas-master/breastcancer_gtex2.top")

# read in TWAS results: TCGA eQTL
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

twas_significant_genes <- twas_significant_genes[order(CHR)]

write.table(twas_significant_genes, "fusion_twas-master/breastcancer_tcga_tumor_ex2.top")

##### FOCUS fine mapping

# extract TWAS significant genes from pos file and create new pos file
twas_results_gtex2 <- fread("fusion_twas-master/breastcancer_gtex2.top")
twas_results_tcga2 <- fread("fusion_twas-master/breastcancer_tcga_tumor_ex2.top")
gtex_pos <- fread("fusion_twas-master/WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue.pos")
tcga_pos <- fread("fusion_twas-master/WEIGHTS/TCGA-BRCA.GE.TUMOR.pos")

w <- which(gtex_pos$ID %in% twas_results_gtex2$ID)
gtex_pos_short <- gtex_pos[w, ]
write.table(gtex_pos_short, file="fusion_twas-master/WEIGHTS/gtex2_twas.pos", row.names=F, col.names=T, sep="\t", quote=F)

w <- which(tcga_pos$ID %in% twas_results_tcga2$ID)
tcga_pos_short <- tcga_pos[w, ]
write.table(tcga_pos_short, file="fusion_twas-master/WEIGHTS/tcga2_twas.pos", row.names=F, col.names=T, sep="\t", quote=F)

# match SNPs for each chromosome
twas_sig_locus_chr <- c(1, 2, 3, 5, 7, 11, 12, 15, 16, 17, 19, 20, 22)
gwas2_header <- colnames(fread("fusion_twas-master/harmonized2.tsv", nrow=0))
gwas2 <- fread("fusion_twas-master/harmonized2.tsv", header=F, col.names=gwas2_header)
for (i in twas_sig_locus_chr) {
  gwas2_i <- copy(gwas2)
  snps <- fread(sprintf("fusion_twas-master/LDREF/1000G.EUR.%s.bim", i), header=F)
  w2 <- which(gwas2_i$SNP %in% snps$V2)
  gwas2_i <- gwas2_i[w2, ]
  gwas2_i$CHR <- i
  gwas2_i$BP <- snps$V4[match(gwas2_i$SNP, snps$V2)]

  write.table(gwas2_i, file=sprintf("harmonized2.sumstats.chr%s", i), row.names=F, col.names=T, sep="\t", quote=F)
}

# run focus finemap (see final.sh)

##### OPTIONAL: RUN JOINT AND CONDITIONAL TESTS (see final.sh)