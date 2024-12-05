#!/bin/bash

for chr in {1..22};
do
	Rscript FUSION.assoc_test.R \
	--sumstats 34737426-GCST90043925-EFO_1000326.rename.h.tsv \
	--weights ./WEIGHTS/TCGA/TCGA-BRCA.TUMOR.pos \
	--weights_dir ./WEIGHTS/TCGA \
	--ref_ld_chr ./LDREF/1000G.EUR. \
	--chr $chr \
	--out TCGA_h_chr${chr}.dat
done
