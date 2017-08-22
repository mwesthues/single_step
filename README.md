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
	* [Preparation](#preparation)
		* [Previous work for this manuscript](#previous-work-for-this-manuscript)
		* [Nested subsampling scheme](#nested-subsampling-scheme)
		* [Prediction template](#prediction-template)
	* [Execution](#execution)
	* [Visualization](#visualization)

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
set up by the Yan lab and are stored under `data/input/maizego`.

The Hohenheim data were generated as part of the publication
[Omics-based Hybrid Prediction in Maize](https://link.springer.com/article/10.1007%2Fs00122-017-2934-0) by Westhues et al. (2017) and can be
found under `data/processed`.


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

## Preparation
### Previous work for this manuscript
In the first draft we wanted to evaluate the influence of the genetic
constitutation of the set of genotypes that has only data on one out of two
predictors by generating core samples  (https://github.com/mwesthues/single_step/commit/8f8b39198be7d5091954f84ec0c5834afd1c3cfa).
We defined core samples as a subset of genotypes that is covered by two
predictors whereas all other genotypes (*i.e.* the complement) are only covered
by a single predictor.
The first application involved assembling core sets of varying sizes by
maximizing the average genetic distance among core set members using the
Modifed Rogers (MR) distance as the criterion ([Thachuk et al. (2009)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-243)).
The size of the core set was varied in increments of ten percentage points and
ranged from 10% to 90% of all inbred lines.
The issue with this approch is that we could not preclude effects of population
structure in our material on predictive abilities.
Therefore, we decided to replace the core sampling procedure by a nested random
subsampling procedure.

### Nested subsampling scheme
At first, a [population structure analysis](analysis/maizego_snp_analyses.R)
was run on the inbred lines and only those that belonged to the clusters `1`
and `4`, which were rather homogeneous, were assigned to scenario `A`.
We defined a second scenario `B` which simply comprised all available inbred
lines and ran it to have predictive abilities for all genotypes in case
population structure was no issue for this material.
Click [here](reports/select_subpopulation.Rmd) for a more detailed report.

To ensure that a potential bias due to population structure would be
represented by an increased standard error of the predictive abilities we
conceived a [nested resampling scheme](analysis/prepare_subsamples.R).
A first level of sampling was applied to each combination of `Material` (*i.e.*
"Hybrid" or "Inbred"), `Scenario` ("A" and "B" for inbred lines; "None" for
hybrids), `Extent` ("Core" if every predictor was available for all genotypes
and "Full for the entire set of genotypes") and `Core_Fraction`, representing
the share of core-set genotypes for which data on all predictors were assumed
to be available.
The notion of `Core_Fraction` may sound contradictory to the concept of "Core"
and "Full" sets of genotypes but note that for the variable `Core_Fraction` we
artificially removed information on a predictor for some genotypes to implore
the influence of genetic space-coverage through incomplete predictors on
predictive ability.
This random subsampling procedure was repeated 20 times and stored as
`Rnd_Level2`.

For all core-set inbred lines we applied a second, nested randomization scheme
where `n * frac` genotypes were declared to have both data on the complete and
the incomplete predictor whereas `n - n * frac` genotypes had only information
on the complete predictor.
Here `frac` denotes the fraction of all genotypes having information on both
predictors.
This second randomization was nested within each randomization stored in
`Rnd_Level2` and applied 20 times per randomization in `Rnd_Level2`.
Therefore, we ended up with 20 sets of predictions for each combination of
`"Material" * "Extent" * "Scenario"` when `frac = 1` and `20 ** 2 = 400`
sets of predictions for each combination of `"Material" * "Extent" *
"Scenario" * "Core_Fraction"` when `frac != 1`.

### Prediction template
In order to use the resources on the server as efficiently as possible,
similar combinations of the various parameters considered (`Trait`,
`Predictor`, `Material`, `Extent`, `Scenario`, `Rnd_Level1`, `Rnd_Level2`,
`Core_Fraction`) were clustered together (coded as `Interval`) in a
[prediction template](analysis/prediction_template.R).

After running some speed tests for some of these combinations a
[script](analysis/automate_moab.R) for the automatic generation of `msub`
commands was built.

## Execution
All [predictions](analysis/prediction.R) were run on the [bwunicluster](https://www.bwhpc-c5.de/wiki/index.php/Category:BwUniCluster)
by entering the [msub commands](analysis/moab_commands.txt) into the terminal
on the cluster.
For information on how to use the `bwunicluster`, in case you have access to
it, can be found [here](https://mwesthues.github.io/bwunicluster.html).



## Visualization


