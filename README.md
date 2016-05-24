This project entails the imputation of gene expression data, which in turn will
be used for the prediction of agronomic traits.
SNP-marker data are available for all inbred lines and will be used as
independent variables throughout the imputation of mRNAs.
Through this intermediary step we should arrive at a genomic prediction
weighted by gene expression data.


## Analysis
1.   [./analysis/snp_qc_pheno.R](./analysis/snp_qc_pheno.R)

SNP quality checks.

2.   [./analysis/snp_impute_pheno.R](./analysis/snp_impute_pheno.R)

SNP imputation using Beagle.

3.   [./analysis/snp_maf_coded.R](./analysis/snp_maf_coded.R)

SNP-coding based on minor-allele frequencies.

4.   [./analysis/equidistant_snps.R](./analysis/equidistant_snps.R)

Select equi-distant SNPs with approximately 10 SNPs per Mbp.

5.   [./analysis/mrna_imputation.R](./analysis/mrna_imputation.R)

Impute mRNAs for inbred lines without mRNA-data using BGLR and genomic kinshipp
information.
