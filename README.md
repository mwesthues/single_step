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

[maizego.R](analysis/maizego.R)

### Agronomic
### Genomic (incl. imputation)
### Gene expression


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

1.   Prepare agronomic 
1.   [common_genotypes](analysis/common_genotypes.R)
2.   [cv_sampling](analysis/cv_sampling.R/)
