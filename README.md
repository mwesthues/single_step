# Table of Contents
<!-- vim-markdown-toc GFM -->
* [Overview](#overview)
	* [Data Preparation](#data-preparation)
	* [UHOH (Hybrids)](#uhoh-hybrids)
		* [Agronomic](#agronomic)
		* [Genomic (incl. imputation)](#genomic-incl-imputation)
		* [Gene expression](#gene-expression)
	* [Yan (Inbreds)](#yan-inbreds)
		* [Agronomic](#agronomic-1)
		* [Genomic](#genomic)
			* [Quality filtering and imputation.](#quality-filtering-and-imputation)
			* [Core set sampling](#core-set-sampling)
			* [Population Structure](#population-structure)
* [Predictions](#predictions)

<!-- vim-markdown-toc -->



# Overview
The purpose of this project is to explore the utility of single-step prediction
of phenotype performance in a set of maize inbred lines and maize hybrids,
respectively.
The methods used herein are based on the seminal papers by [Legarra et al. (2009)](http://www.sciencedirect.com/science/article/pii/S0022030209707933)
and [Christensen and Lund (2010)](https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-42-2), with extensions introduced by
[Fernando, Dekkers and Garrick (2014)](https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-46-50).






## Data Preparation
We analyze hybrid maize data from the public breeding program of the University
of Hohenheim as well as [data on inbred lines](http://www.maizego.org/Resources.html) from the lab of doctor Jianbing Yan.

The data on inbred lines were downloaded on 2017-01-31 from the [Baidu Cloud](https://pan.baidu.com/s/1eQH3hfW#list/path=%2F)
set up by the Yan lab and are stored under [data/input/maizego](data/input/maizego).

The Hohenheim data were generated as part of the publication
[Omics-based Hybrid Prediction in Maize](https://link.springer.com/article/10.1007%2Fs00122-017-2934-0) by Westhues et al. (2017) and can be
found under [data/processed](data/processed/).


Run all scripts strictly in the following order in order to circumvent any
issues with dependencies:

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

