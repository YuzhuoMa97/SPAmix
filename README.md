# SPAmix
A scalable, accurate, and universal analysis framework using individual-level allele frequency for large-scale genetic association studies in an admixed population.


## Table of contents
  * [Software dependencies and operating systems](#Software-dependencies-and-operating-systems)
  * [How to install and load package](#How-to-install-and-load-package-for-SPAmix-analysis)
  * [Introduction of SPAmix](#Introduction-of-SPAmix)
  * [Summary of key features for our proposed scalable methods](#Summary-of-key-features-for-our-proposed-scalable-methods)
  * [Reproducibility](#Reproducibility)
      * [UK Biobank data analysis](#UK-Biobank-data-analysis)
      * [Simulation studies](#Simulation-studies)
  * [Reference](#Reference)


## Software dependencies and operating systems

- SPAmix has been implemented in GRAB package. The package has been thoroughly examined and validated on both Linux and Windows operating systems.

- Currently, this R package supports three formats for genotype input: the R matrix (Rdata) format, the [**PLINK**](https://www.cog-genomics.org/plink2/) format, and the [**BGEN**](https://www.chg.ox.ac.uk/~gav/bgen_format/index.html) format.

## How to install and load this package for SPAmix analysis
```
# ========================================================================
# SPAmix Analysis Pipeline using GRAB Package
# ========================================================================

# -------------------------------------------
# Package Installation (Uncomment if needed)
# -------------------------------------------
# For Linux/macOS systems:
# library(remotes)  
# install_github("GeneticAnalysisinBiobanks/GRAB")  

# For Windows systems (requires Rtools):
# library(remotes)
# options(download.file.method = "wininet")
# install_github("GeneticAnalysisinBiobanks/GRAB", INSTALL_opts=c("--no-multiarch"))

# Load core package
library(GRAB)
```

- Please refer to the [GRAB documentation](https://wenjianbi.github.io/grab.github.io/docs/approach_SPACox.html) for more details on how to conduct GWAS using this package.

- Please do not hesitate to contact yuzhuoma@stu.pku.edu.cn or wenjianb@pku.edu.cn if you meet any problem. Suggestions or comments are also welcome.

## Introduction of SPAmix 

### SPAmix is designed for conducting GWASs in a study cohort including multiple subpopulations and/or admixed populations. 

**SPAmix is applicable to a wide range of complex traits with intricate structures, including time-to-event, ordinal categorical, binary, quantitative, longitudinal, and other complex traits.** The framework involves two main steps:

- Step 1: SPAmix requires fitting a null model to calculate model residuals, in which confounding factors such as age, sex, SNP-derived principal components (PCs), and other confounders are incorporated as covariates. The null model specification and the corresponding model residuals can vary depending on the type of trait. In the online Methods section, we demonstrated regression models to fit binary, quantitative, time-to-event, ordinal, and longitudinal traits, and the corresponding model residuals. 

- Step 2: SPAmix approach associates the traits of interest to single genetic variant by approximating the null distribution of score statistics. To characterize the diversity of allele frequencies (AFs) of this genetic variant for different ancestries, we assume that individuals come from different populations with varying AF. Instead, linear and logistic regressions were used to estimate individual-specific AFs, in which raw genotypes are outcomes and SNP-derived principal components are model predictors. SPAmix uses a hybrid strategy including both normal distribution approximation and SPA to approximate the distribution of score statistics under the null hypothesis. In addition, model residuals are categorized to outliers and non-outliers based on interquartile range, and a partial normal distribution is used to reduce computational burden. SPAmix supports simultaneously analyzing multiple traits whose types are the same or different, which can avoid repeated calculation of the AFs estimation. 



# Introduction of Retrospective saddlepoint approximation (Retrospective SPA) approach in GWAS

Retrospective saddlepoint approximation (Retrospective SPA) is a method applied in GWAS. For a score statistic (S=G^TR), retrospective SPA strategy considers genotypes as random variables and approximate the distribution of score statistics conditional on phenotype and covariates.

Suggestions or comments on retrospective saddlepoint approximation methods are also welcome.

# SPAmix can control for population admixture

In the SPAmix paper (**SPAmix: A scalable, accurate, and universal analysis framework for large-scale genetic association studies in admixed populations**), we proposed SPAmix framework that is applicable to admixed populations and can incorporate local ancestry information into analyses. 

