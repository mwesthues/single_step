# Table of Contents
<!-- vim-markdown-toc GFM -->
* [Overview](#overview)
* [Data Preparation](#data-preparation)
	* [Agronomic Data](#agronomic-data)
	* [Genomic Data](#genomic-data)
	* [Transcriptomic Data](#transcriptomic-data)
* [Predictions](#predictions)
	* [Preparation](#preparation)
		* [Prediction template](#prediction-template)
	* [Execution](#execution)
	* [Bootstrap](#bootstrap)
	* [Visualization](#visualization)

<!-- vim-markdown-toc -->



# Overview
The purpose of this project is to explore the utility of single-step prediction
of phenotype performance in a set of maize inbred lines and maize hybrids,
respectively.
The methods used herein are based on the seminal papers by [Legarra et al. (2009)](http://www.sciencedirect.com/science/article/pii/S0022030209707933)
and [Christensen and Lund (2010)](https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-42-2), with extensions introduced by
[Fernando, Dekkers and Garrick (2014)](https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-46-50).






# Data Preparation
We analyze hybrid maize data from the public breeding program of the University
of Hohenheim as well as [data on inbred lines](http://www.maizego.org/Resources.html) from the lab of doctor Jianbing Yan.

The data on inbred lines were downloaded on 2017-01-31 from the [Baidu Cloud](https://pan.baidu.com/s/1eQH3hfW#list/path=%2F)
set up by the Yan lab. Specifically, the following files were downloaded and
stored under `data/input/maizego`:

*   The file `Expression_final_QQNormalized.tar.gz` from the [Baidu
    Cloud](https://pan.baidu.com/s/1eQH3hfW#list/path=%2FMaizeGo%20Resources%2FTranscriptomic&parentPath=%2F).
*   The file `blup_traits_final.csv` from the [Baidu
    Cloud](https://pan.baidu.com/s/1eQH3hfW#list/path=%2FMaizeGo%20Resources%2FPhenotype%2FAgronomic&parentPath=%2F).
*   The file `513allchr-imputation.zip` from the [Baidu
    Cloud](https://pan.baidu.com/s/1eQH3hfW#list/path=%2FMaizeGo%20Resources%2FGenotype%2FSNPs_513lines_0.55M_MAF0.05&parentPath=%2F).
*   The file `513lines_27229snps_kinship_110608.txt` from the [Baidu
    Cloud](https://pan.baidu.com/s/1eQH3hfW#list/path=%2FMaizeGo%20Resources%2FGenotype&parentPath=%2F).

The Hohenheim data were generated as part of the publication
[Omics-based Hybrid Prediction in Maize](https://link.springer.com/article/10.1007%2Fs00122-017-2934-0) by Westhues et al. (2017) and can be
found [here](data/processed/).


Run all scripts strictly in the following order in order to circumvent any
issues with dependencies:

## Agronomic Data
For the maize inbred lines we focused on the following traits that were
previously described in [Guo et al. (2016)](https://link.springer.com/article/10.1007/s00122-016-2780-5):

*   `Cob weight` (CW)
*   `Days to silking` (DS)
*   `Ear diameter` (ED)
*   `100-grain weight` (GW)
*   `Kernel width` (KW)
*   `Plant height` (PH)

In the case of maize hybrids we analyzed the traits

*   `Dry matter yield` (DMY)
*   `Dry matter content` (DMC)
*   `Fat` (FAT)
*   `Protein` (PRO)
*   `Starch` (STA)
*   `Sugar` (SUG)
*    Lignin (ADL)

The preparation of the phenotypic data for maize hybrids can be reproduced
[here](analysis/uhoh_data_preparation.R).


## Genomic Data
The preparation of genomic marker data for [hybrids](analysis/uhoh_snp_imputation.R)
and [inbred lines](analysis/maizego_snp_imputation.R), respectively, involved
the following five steps (for hybrid data these were applied separately to each
of the two heterotic groups):

1.   Remove loci with a callfrequency $\geq 0.95$.

2.   Remove loci with heterozygosity $\geq 0.05$.

3.   Remove loci with minor allele frequency $\geq 0.05$.

4.   Impute remaining missing values using BEAGLE using 25 iterations and 20
     samples per iterations.

5.   Remove copies (i.e. loci in perfect LD with a previous one) of marker
     loci.


## Transcriptomic Data
For both material types transcriptomic data were already sufficiently
pre-processed.
Their assembly, with respect to the set of genotypes used in this study, can
be found under the two following links
([hybrids](analysis/uhoh_data_preparation.R), [inbreds](analysis/maizego_gene_expression.R)).


# Predictions
## Preparation
Test and training sets and their corresponding predictor data matrices were
generated in a [separate script](analysis/prepare_subsamples.R) prior to
the predictions.

### Prediction template
In order to use the resources on the server as efficiently as possible,
similar combinations of the various parameters considered (`Trait`,
`Predictor`, `Material`, `Extent`, `Scenario`, `Rnd_Level1`, `Core_Fraction`)
were clustered together (coded as `Interval`) in a
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

After running the prediction, the results from all clusters were [concatenated](analysis/cat_predictions.R).


## Bootstrap
Separately for each combination of `Material`, `Extent`, `Scenario`, `Trait`,
`Core_Fraction`, `Predictor` and `Rnd_Level1` a [bootstrap](analysis/bootstrap_predictions.R)
of 10,000 runs was applied to a data frame comprising observed and predicted
phenotypic values.
This yielded 10,000 bootstrap predictive abilities for each of the
aforementioned combinations.


## Visualization
Predictive abilities were plotted using [this script](analysis/visualize_predictions.R).

