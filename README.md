This project entails the imputation of gene expression data, which in turn will
be used for the prediction of agronomic traits.
SNP-marker data are available for all inbred lines and will be used as
independent variables throughout the imputation of mRNAs.
Through this intermediary step we should arrive at a genomic prediction
weighted by gene expression data.


## Analysis
### Common genotypes
[./analysis/common_genotypes.R](./analysis/common_genotypes.R)

Determine for which genotypes (inbred lines as well as hybrids) all data are
available.

### Genomic data
[./analysis/snp_preparation.R](./analysis/snp_preparation.R)

Apply genomic data pre-processing such as the imputation of missing values.

### Hybrid prediction
[./analysis/prediction.R](./analysis/prediction.R)

Run hybrid prediction models with the `BGLR` package using either leave-one-out
cross validation (LOOCV) or customized CV schemes from 
[./analysis/cv_sampling.R].


### Model evaluation
[./analysis/pheno_predability.R](./analysis/pheno_predability.R)

Compute and plot predictive abilities.




## Helper scripts
### Predictive ability
[./analysis/helper_predability.R](./analysis/helper_predability.R)

Determine predictive abilities for single-step BLUP-based predictions.
Predictive abilities are plotted separately for hybrids for which mRNA values
were imputed for none, one or both parents, each for CV1000 and LOOCV.

### CV1000 composition
[./analysis/helper_loocv_trn_size.R](./analysis/helper_loocv_trn_size.R)

Determine how similar the individual cross-validations rounds in CV1000 are.




## Additional scripts (not required)
### CV scheme generation
[./analysis/cv_sampling.R](./analysis/cv_sampling.R)

Generate custom cross-validation schemes as outlined in Westhues et al. (2016).



