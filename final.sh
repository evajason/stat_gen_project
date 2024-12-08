# run FUSION on all chromosomes

# first GCST sumstats & GTEx eQTL
for xx in {1..22}
do
Rscript FUSION.assoc_test.R \
	--sumstats GCST90043925_harmonized.tsv.gz \
	--weights ./WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue.pos \
	--weights_dir ./WEIGHTS/ \
	--ref_ld_chr ./LDREF/1000G.EUR. \
	--chr ${xx} \
	--out breastcancer_chr${xx}_gtex.dat
done

# second GCST sumstats & GTEx eQTL
for xx in {1..22}
do
Rscript FUSION.assoc_test.R \
    --sumstats GCST90454347.v2.tsv.gz \
    --weights ./WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue.pos \
    --weights_dir ./WEIGHTS/ \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr ${xx} \
    --out breastcancer_chr${xx}_gtex2.dat
done

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
# first GCST sumstats & TCGA eQTL
for xx in {1..22}
do
Rscript FUSION.assoc_test.R \
    --sumstats GCST90043925_harmonized.tsv.gz \
    --weights ./WEIGHTS/TCGA-BRCA.GE.TUMOR.pos \
    --weights_dir ./WEIGHTS/ \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr ${xx} \
    --out breastcancer_chr${xx}_tcga_tumor_exp.dat
done

# second GCST sumstats & TCGA eQTL
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

# including coloc
Rscript FUSION.assoc_test.R \
	--sumstats GCST90043925_harmonized.tsv.gz \
	--weights ./WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue.pos \
	--weights_dir ./WEIGHTS/ \
	--ref_ld_chr ./LDREF/1000G.EUR. \
	--chr 2 \
	--out breastcancer_chr2_gtex.dat \
	--coloc_P 0.05 

# run joint and conditional tests
# how much GWAS signal remains after the association of the functional is removed?
Rscript FUSION.post_process.R \
	--sumstats GCST90043925_harmonized.tsv.gz \
	--input breastcancer_gtex.top \
	--out breastcancer_chr7_gtex.top.analysis \
	--ref_ld_chr ./LDREF/1000G.EUR. \
	--chr 7 \
	--plot --locus_win 100000

# FOCUS fine mapping
focus import fusion_twas-master/WEIGHTS/GTExv8.EUR.Breast_Mammary_Tissue.pos fusion --tissue breastmammary --name GTEx --assay rnaseq --output fusion_breastmammary

focus import fusion_twas-master/WEIGHTS/TCGA-BRCA.GE.TUMOR.pos fusion --tissue breastmammary --name TCGA --assay rnaseq --output fusion_tcga_breastmammary

focus finemap harmonized.sumstats.chr2 fusion_twas-master/LDREF/1000G.EUR.2 fusion_breastmammary.db \
	--chr 2 \
	--prior-prob "gencode38" \
	--locations 38:EUR \
	--out ./focus_result_bc2

for i in 2 7 11
do
focus finemap harmonized.sumstats.chr${i} fusion_twas-master/LDREF/1000G.EUR.${i} fusion_breastmammary.db \
	--chr ${i} \
	--prior-prob "gencode38" \
	--locations 38:EUR \
	--out ./focus_result_bc${i}
done

for i in 5 11 22
do
focus finemap harmonized.sumstats.chr${i} fusion_twas-master/LDREF/1000G.EUR.${i} fusion_tcga_breastmammary.db \
	--chr ${i} \
	--prior-prob "gencode38" \
	--locations 38:EUR \
	--out ./focus_result_bc_tcga${i}
done
