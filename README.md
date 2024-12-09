## Introduction

This is a final project for Statistical Genomics (DSC 291) taught at UCSD Fall 2024

Given the impact of genetics on breast cancer, we would like to see if TWAS and GWAS have any shared loci of expressed genes and variants associated with breast cancer. Going beyond associations, we would like to see if there are any “causal” genes identified using fine-mapping.

## Set Up:

1. Install FUSION: Follow the installation instructions according to the guidelines on their [repository](http://gusevlab.org/projects/fusion/).
2. Instal FOCUS: Follow the installation instructions according to the guidelines on their [repository](https://github.com/mancusolab/ma-focus).

3. Download the necessary files (remember to cd in to the `fusion_twas-master` directory):

- [Gene Expression Weights (GTEx)](https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Breast_Mammary_Tissue.tar.gz)
- [Gene Expression Weights (TCGA breast tumor expression)](http://gusevlab.org/projects/fusion/weights/GusevLawrenson_2019_NG/TCGA-BRCA.GE.TUMOR.tar.bz2)
- Harmonized GWAS summary statistics for breast cancer from the GWAS Catalog created by EMBL-EBI, study [GCST90454347](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454347/harmonised/GCST90454347.h.tsv.gz)

4. Using the instructions on the [FUSION repository](http://gusevlab.org/projects/fusion/), set up the `WEIGHTS/` and `LDREF/` directories.

## Run the analysis:

3. Modify the GWAS sumstats files for FUSION

- Rename the columns so that we have `SNP`, `A1`, and `A2`, and calculate the `Z` column.

  Use the make_z_col.R script in this repository.

`Rscript make_z_col.R GCST90454347.h.tsv.gz` 

The script will take a GWAS summary stats file and write it out as a tsv with the correct columns.

4. Run FUSION on both predictive models (GTEx and TCGA)

- See `final.sh` for code chunks. Please ensure that your edited sumstats file is named `harmonized2.tsv`. 
- After running FUSION on all chromosomes, see `final.R` to aggregate the TWAS significant genes across all chromosomes. 

5. Run FOCUS (install FOCUS [here](https://github.com/mancusolab/ma-focus))

- See `final.R` for fine-mapping setup. Run focus finemap from `final.sh`.

6. Compare GWAS and TWAS significant results with FOCUS results

- See `final_proj.Rmd` for code chunks.
