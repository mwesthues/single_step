Run all scripts strictly in the following order in order to circumvent any
issues with dependencies:

# Data preparation
## UHOH (Hybrids)
### Agronomic
### Pedigree
### Genomic (incl. imputation)
### Gene expression





## Yan (Inbreds)
Filter agronomic, genomic and transcriptomic data for tropical and subtropical
lines.

[maizego_tst_data_subsetting.R](analysis/maizego_tst_data_subsetting.R)



### Agronomic
1.   Distribution of agronomic data.

2.   Skewness in agronomic data (permutation test).

3.   Correlation among agronomic traits.

[maizego_agronomic_data.R](analysis/maizego_agronomic_data.R)


### Genomic
#### Quality filtering and imputation.

1.   Remove loci with a callfrequency $\geq 0.95$.

2.   Remove loci with heterozygosity $\geq 0.05$.

3.   Remove loci with minor allele frequency $\geq 0.05$.

4.   Impute remaining missing values using BEAGLE using 25 iterations and 20
     samples per iterations.

5.   Remove copies (i.e. loci in perfect LD with a previous one) of marker
     loci.

[maizego_snp_imputation.R](analysis/maizego_snp_imputation.R)

6.   Principal component analysis.



#### Core set sampling



### Gene expression
1.   Scaling, centering, Box-Cox transformation

2.   Principal component analysis

[maizego_mrna_preprocessing.R](analysis/maizego_mrna_preprocessing.R)





Create 'predictor_data/' directory where, separately for 'maizego' and 'uhoh',
predictor data for the predictions are stored.


# Predictions
## UHOH (Hybrids)
## Scenario 1
## Scenario 2
## Scenario 3

## Yan (Inbreds)
## Scenario 1
## Scenario 2
## Scenario 3

