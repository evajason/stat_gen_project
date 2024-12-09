# Set Up:

1. Install FUSION: Follow the installation instructions according to the guidelines on their [repository](http://gusevlab.org/projects/fusion/).

2. Download the necessary files (remember to cd in to the `fusion_twas-master` directory):

- [Gene Expression Weights (GTEx)](https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Breast_Mammary_Tissue.tar.gz)
- [Gene Expression Weights (TCGA breast tumor expression)](http://gusevlab.org/projects/fusion/weights/GusevLawrenson_2019_NG/TCGA-BRCA.GE.TUMOR.tar.bz2)
- Harmonized GWAS summary statistics for breast cancer from the GWAS Catalog created by EMBL-EBI, study [GCST90454347](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454347/harmonised/GCST90454347.h.tsv.gz)

(Following the steps and code in `final.R` and `final.sh`...) 

3. Modify the GWAS sumstats files for FUSION

- Rename the columns so that we have `SNP`, `A1`, and `A2`, and calculate the `Z` column. 

4. Run FUSION on both predictive models (GTEx and TCGA)

- See `final.sh` for code chunks. Please ensure that your edited sumstats file is named `harmonized2.tsv`. 
- After running FUSION on all chromosomes, see `final.R` to aggregate the TWAS significant genes across all chromosomes. 

5. Run FOCUS (install FOCUS [here](https://github.com/mancusolab/ma-focus))

- See `final.R` for fine-mapping setup. Run focus finemap from `final.sh`.

6. Compare GWAS and TWAS significant results with FOCUS results

- See `final_proj.Rmd` for code chunks.
