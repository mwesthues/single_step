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
