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

[maizego_snp_analyses.R](analysis/maizego_snp_analyses.R)


#### Core set sampling
##### Scenario 1
Scenario 1 does not involve any core set sampling because it simply uses all
available SNP data (211 genotypes) and all available transcriptomic data (149
genotypes).
However, it does involve subsampling of the SNP data for the direct comparison
of the predictive ability obtained for 149 genotypes with transciptomic data
and the one obtained for the same set of genotypes using genomic data.


##### Scenario 2
Goal: Quantify the influence of SNP data on the predictive ability of the
combination of genomic with transcriptomic data in single step prediction.

1.   Reduce the SNP data to 149 inbred lines, which are also covered by mRNA 
     data.

2.   Alter the fraction of genotypes covered by SNP data from 100% to 10% in
     increments of 10 percentage points while keeping the number of genotypes
     covered by the transcriptomic data fixed.

[maizego_corehunter.R](analysis/maizego_corehunter.R)


##### Scenario 3
This scenario is based on scenario 2 but it respects the set of genotypes
covered by trancsriptomic data.
This means that, as long as the fraction of genotypes covered by SNP data is at
least as large as the set of genotypes covered by transcriptomic data, all
genotypes covered by transcriptomic data must be included in the core set.
As soon as the set of genotypes covered by SNP data is smaller than the one
covered by mRNA data, we sample the core set only from the set of genotypes
covered by mRNA data.

> I'm still undecided whether this adds anything useful to our study compared
> to scenario 2.

[maizego_corehunter.R](analysis/maizego_corehunter.R)


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

