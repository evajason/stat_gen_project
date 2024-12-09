# run FUSION on all chromosomes

# GTEx eQTL weights
for xx in {1..22}
do
Rscript FUSION.assoc_test.R \
    --sumstats harmonized2.tsv \
    --weights ./WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue.pos \
    --weights_dir ./WEIGHTS/ \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr ${xx} \
    --out breastcancer_chr${xx}_gtex2.dat
done

# TCGA eQTL weights
for xx in {1..22}
do
Rscript FUSION.assoc_test.R \
    --sumstats harmonized2.tsv \
    --weights ./WEIGHTS/TCGA-BRCA.GE.TUMOR.pos \
    --weights_dir ./WEIGHTS/ \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr ${xx} \
    --out breastcancer_chr${xx}_tcga_tumor_exp2.dat
done

# FOCUS fine mapping
focus import fusion_twas-master/WEIGHTS/gtex2_twas.pos fusion --tissue breastmammary --name GTEx --assay rnaseq --output fusion_gtex_breastmammary

focus import fusion_twas-master/WEIGHTS/tcga2_twas.pos fusion --tissue breastmammary --name TCGA --assay rnaseq --output fusion_tcga_breastmammary

# GTEx eQTL
for i in 1 3 5 7 12 15 17 19 22
do
focus finemap harmonized2.sumstats.chr${i} fusion_twas-master/LDREF/1000G.EUR.${i} fusion_gtex_breastmammary.db \
    --chr ${i} \
    --prior-prob "gencode38" \
    --location 38:EUR \
    --out ./focus_results_bc_gtex2_chr${i}
done

# TCGA eQTL
for i in 1 3 5 12 15 16 20 22
do
focus finemap harmonized2.sumstats.chr${i} fusion_twas-master/LDREF/1000G.EUR.${i} fusion_tcga_breastmammary.db \
    --chr ${i} \
    --prior-prob "gencode38" \
    --locations 38:EUR \
    --out ./focus_result_bc_tcga${i}
done

# OPTIONAL: run conditional tests; how much GWAS signal remains after the association of the functional is removed?

# GTEx eQTL FOCUS causal genes
for xx in 3 12 15 
do
Rscript FUSION.post_process.R \
    --sumstats harmonized2.tsv \
    --input breastcancer_gtex2.top \
    --out breastcancer_chr${xx}_gtex2.top.analysis \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr ${xx} \
    --plot --locus_win 100000
done

# TCGA eQTL FOCUS causal genes
for xx in 5 11 22
do
Rscript FUSION.post_process.R \
    --sumstats GCST90043925_harmonized.tsv.gz \
    --input breastcancer_tcga_tumor_ex.top \
    --out breastcancer_chr${xx}_tcga_tumor_ex.top.analysis \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr ${xx} \
    --plot --locus_win 100000
done
