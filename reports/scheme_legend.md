# a
Level 1 sampling scheme.
During the first round of sampling, 70% of the genotypes were drawn at random
(black symbols) and assigned to the training set (TRN) whereas phenotypes for
the other 30% of genotypes (yellow symbols) were set to missing for this
particular random sample.
The random samples drawn at this level were subsequently denoted as s~1~ to
s~20~.

# b
Level 2 sampling scheme.
For maize inbred lines from the diversity panel, some genotypes had both
genotypic and transcriptomic predictor data (black symbols).
For each level of s from the first round of random sampling, a second round
of random sub-sampling was applied where x% of the genotypes from the training
set of s were declared to have both, genotypic and transcriptomic data (black
symbols) whereas 1-x% of the genotypes were declared to have only genotypic
information (gray symbols).
This nested subsampling process was repeated 20 times and for nine different
values of x (0.1, 0.2, 0.3, ..., 0.9).

# c
Principal component analysis of maize inbred lines from the diversity panel.
A STRUCTURE analysis with K=4 putative ancestral populations was run to further
explore population structure.
Genotypes were colored based on their highest share of membership with any of
the four putataive ancestral populations.

# d
Summary of both levels of random sampling.
The level 1 randomization refers to **a** and emphasizes that the first level of
randomization was nested within scenarios but identical across all traits.
Scenario A entails only genotypes from clusters 1 and 4 in **c**.
Scenario B and C comprise all genotypes in **c** where scenario C additionally
accounts for the four different clusters by modeling them as fixed effects in
the predictions.
Level 2 randomization refers to **b** and emphasizes that, for each scenario,
each set from level 1 and for each x-value, 20 different subsets were sampled.


# e
Leave-one-out cross-validation (LOOCV) scheme for maize diversity panel inbred
lines.
For each of the subsets from level 1 (and level 2) LOOCV was applied.
Thereto, n TST were generated (yellow symbols) where all genotypes eligible for
TRN membership based on **a** were assigned to TRN (black symbols).

# f
LOOCV scheme for maize hybrids.
This scheme is based solely on level 1 randomization in **a** and, for each
level of s, places only T0 hybrids, which do not share a common parent with a
TST hybrid, inside the TRN.

