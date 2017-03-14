Run all scripts strictly in the following order in order to circumvent any
issues with dependencies:

# Data preparation
## UHOH (Hybrids)
### Agronomic
[agronomic_data.Rmd](reports/agronomic_data.Rmd)

### Genomic (incl. imputation)

1.   Remove loci with a callfrequency $\geq 0.95$.

2.   Remove loci with heterozygosity $\geq 0.05$.

3.   Remove loci with minor allele frequency $\geq 0.05$.

4.   Impute remaining missing values using BEAGLE using 25 iterations and 20
     samples per iterations.

5.   Remove copies (i.e. loci in perfect LD with a previous one) of marker
     loci.

[snp_preparation.R](analysis/snp_preparation.R)


### Gene expression

[transcriptomic_data.Rmd](reports/transcriptomic_data.Rmd)




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




#### Core set sampling
Goal: Quantify the influence of SNP data on the predictive ability of the
combination of genomic with transcriptomic data in single step prediction.

1.   Reduce the SNP data to 149 inbred lines, which are also covered by mRNA 
     data.

2.   Alter the fraction of genotypes covered by SNP data from 100% to 10% in
     increments of 10 percentage points while keeping the number of genotypes
     covered by the transcriptomic data fixed.


#### Population Structure

1.   Principal component analysis highlighting genotypes with mRNA data only
     versus genotypes with mRNA and genomic information.

2.   Admixture analysis.

[maizego_snp_analyses.R](analysis/maizego_snp_analyses.R)


3.  PCA highlighting membership with a core set.

[maizego_corehunter.R](analysis/maizego_corehunter.R)






# Predictions

[prediction.R](analysis/prediction.R)

