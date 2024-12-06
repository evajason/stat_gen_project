#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# 34737426-GCST90043925-EFO_1000326.h.tsv.gz
# 34737426-GCST90043925-EFO_1000326.rename.h.tsv

# This script takes a GWAS file with columns called 'beta' and 'standard error'
# and out puts the a file with these columns renamed and a column called Z for 
# the z-score
library(dplyr)
library(data.table)

df = fread(args[1]) %>%
  mutate(Z = beta/standard_error) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele)
write.table(df, args[2])