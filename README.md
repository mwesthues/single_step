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

5.   [./analysis/legarra_kernels.R](./legarra_kernels.R)

Generate kernels according to equation (4) in Legarra et al. (2009)

6.   [./analysis/CV1000_Sampling.R](./analysis/CV1000_Sampling.R)

Generate a CV1000 (1,000 rounds of cross validation) scheme.

7.   [./analysis/SS-BLUP.R](./analysis/SS-BLUP.R)

Run single step regression BLUP.



#### Previous scripts
a.   [./analysis/pheno_prediction.R](./analysis/pheno_prediction.R)

Prediction of agronomic traits in maize hybrids based on mRNAs, which were
imputed for genotypes for which only genomic data were available.

b.   [./analysis/pheno_predability.R](./analysis/pheno_predability.R)

Compute the predictive abilities from the cross-validation runs and plot them.
