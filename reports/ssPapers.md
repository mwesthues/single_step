---
title: "single-step BLUP Paper"
output: 
  html_document:
    toc: true
    toc_float: true
    smooth_scroll: true
    toc_depth: 3
    number_sections: true
    theme: lumen
    highlight: tango
    fig_width: 8
    fig_height: 5
---



# Legarra-2009-JournalDairyScience
## Goal
It is proposed to condition the genetic value of ungenotyped animals on the
genetic value of genotypes animals via the selection index (e.g. pedigree
information), and then use the genomic relationship matrix for the latter. 
This results in a joint distribution of genotypes and ungetnoyped genetic
values with a pedigree-genomic relationship matrix $\mathbf{H}$.
In this matrix, genomic information is trnsmitted to the covariances among all
ungenotypes individuals.

## Introduction
Although these methods [whole genome regression] are very promising for animal
breeding, genotyping is not feasible for an entire population because of its
high cost or logistical constraints.


## Methods
### Covariance matrix of breeding values including genomic information
VanRaden [@VanRaden2008] discussed how the expectation of $\mathbf{G}$ is
$\mathbf{A}$, the regular numerator relationship matrix, and that $\mathbf{G}$
represents observed, rather than average, relationships.
Therefore it accounts for Mendelian samplings (*i.e.* it can distinguish
full-sibs and unknown or far relationships.

#### Plug-in G
A simple way to use $\mathbf{G}$ is to plug it into $\mathbf{A}$; this
results in the following modified $\mathbf{A}$:

$$
  A_{g} = \begin{bmatrix}
            \mathbf{A}_{11} & \mathbf{A}_{12} & \mathbf{A}_{13} \\
            \mathbf{A}_{21} & G & \mathbf{A}_{23} \\
            \mathbf{A}_{31} & \mathbf{A}_{32} & \mathbf{A}_{33} \\
          \end{bmatrix}
$$
Matrix $\mathbf{A}_{g}$ is simple to use but not properly constructed.
The use of $\mathbf{G}$ potentially modified covariances in ancestors and
descendants of genotypes animals.
For example, assume two full-sibs in the genotypes animals whose genomic
relationsihp is 0.6.
By using $mathbf{A}_{g}$, it is assumed that aaaverage relationship among their
daughters is 0.25, whreas in fact it is 0.3."

It can be verified by small numerical examples that $\mathbf{A}_{g} is indefinite
(*i.e.* some eigenvalues are negative and some positive). [...]
A correct statistical inference can only be made if the covariance matrix is
positive or semi-positive definite.



# Fernando-2014-GSE
## Introduction
When $n_{2}$ (the number of genotypes animals) is a few thousands, SSBV-BLUP 
provides an elegant and convenient method to estimate breeding values (BV) that 
combines the available phenotype, pedigree and SNP data.
Due to widespread adoption of gneotpying in livestock, however, $n_{2}$ is 
becoming too large for SSBV-BLUP to remain computationally feasible much
longer.
As discussed by Stranden and Garrick (2009), when the number of gneotypes
individuals exceeds the number of marker covariates, use of marker effect
models (MEM), wich do not require computing $\mathbf{G}$ or its inverse, will lead 
to more efficient calculations.
Our extended MEM will enable BLUP evaluations without having to compute
$\mathbf{G}$ or its inverse, while combining information from genotyped and
non-gneotyped animals.

## Methods
##### Marker effect models
General model:

$$
  \mathbf{y} = \mathbf{X \beta} + \mathbf{T \alpha} + \mathbf{e},
$$ {#eg:General_Model}
where $\mathbf{y}$ is the vector of trait phenotypes, \mathbf{X} is a known
incidence matrix that relates the vector of non-genetic, "fixed" effects to
$\mathbf{y}$, $\mathbf{T} = \mathbf{M} - E(\mathbf{M})$, $\mathbf{M}$ is a
matrix of marker covariates, $\mathbf{\alpha}$ is the vector of random, partial
regression coefficients of the marker covariates, and $\mathbf{e}$ is a vector
of residuals.

### Breeding value models
A model that is equivalent to ([@eq:General_Model]) can be written as:
$$
  \mathbf{y} = \mathbf{X \beta} + \mathbf{g} + \mathbf{e},
$$ {#eg:BV_Model}
where $\mathbf{g} = \mathbf{T \alpha}$, and $\mathbf{T} = \mathbf{M} -
E(\mathbf{M})$ is the matrix of centered marker covariates.


##### Theory underlying SSBV-BLUP

