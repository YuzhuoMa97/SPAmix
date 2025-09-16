# SPAmix
A scalable, accurate, and universal analysis framework using individual-level allele frequency for large-scale genetic association studies in an admixed population.


## Table of contents
  * [Software dependencies and operating systems](#Software-dependencies-and-operating-systems)
  * [How to install and load package](#How-to-install-and-load-package)
  * [Introduction of SPAmix](#Introduction-of-SPAmix)
  * [Summary of key features for our proposed scalable methods](#Summary-of-key-features-for-our-proposed-scalable-methods)
  * [Reproducibility](#Reproducibility)
      * [UK Biobank data analysis](#UK-Biobank-data-analysis)
      * [Simulation studies](#Simulation-studies)
  * [Reference](#Reference)


## Software dependencies and operating systems

- The package has been thoroughly examined and validated on both Linux and Windows operating systems.

- Currently, this R package supports three formats for genotype input: the R matrix (Rdata) format, the [**PLINK**](https://www.cog-genomics.org/plink2/) format, and the [**BGEN**](https://www.chg.ox.ac.uk/~gav/bgen_format/index.html) format.

- In the near future, this R package is planned to be rewritten using Rcpp code to improve its performance and efficiency.

## How to install and load this package
```
library(devtools)  # author version: 2.4.5
install_github("YuzhuoMa97/SPAGxECCT")
library(SPAGxECCT)
?SPAGxECCT  # manual of SPAGxECCT package
```
- Current version is 1.1.0. For older version and version update information, plesase refer to OldVersions.  

- Please refer to the [SPAGxE documentation](https://yuzhuoma97.github.io/RetroSPAgwas.github.io/docs/approach_GxE.html) for more details on how to conduct GxE analysis using this package.

- Please do not hesitate to contact yuzhuoma@stu.pku.edu.cn or wenjianb@pku.edu.cn if you meet any problem. Suggestions or comments are also welcome.

## Introduction of SPAmix 

### SPAGxE<sub>CCT</sub> is a scalable and accurate G×E analytical framework that accounts for unbalanced phenotypic distribution

**SPAGxE<sub>CCT</sub> is applicable to a wide range of complex traits with intricate structures, including time-to-event, ordinal categorical, binary, quantitative, longitudinal, and other complex traits.** The framework involves two main steps:

- Step 1: SPAGxE<sub>CCT</sub> fits a covariates-only model to calculate model residuals. These covariates include, but are not limited to, confounding factors such as age, sex, SNP-derived principal components (PCs), and environmental factors. The specifics of the model and residuals vary depending on the trait type. Since the covariates-only model is genotype-independent, it only needs to be fitted once across a genome-wide analysis.

- Step 2: SPAGxE<sub>CCT</sub> identifies genetic variants with marginal G×E effects on the trait of interest. First, marginal genetic effects are tested using score statistics. If the marginal genetic effect is not significant, S<sub>G×E</sub> is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E</sub> is updated to genotype-adjusted test statistics. To balance computational efficiency and accuracy, SPAGxE<sub>CCT</sub> employs a hybrid strategy combining normal distribution approximation and saddlepoint approximation (SPA) to calculate p-values, as used in previous studies such as [SAIGE](https://saigegit.github.io/SAIGE-doc/) and [SPAGE](https://github.com/WenjianBI/SPAGE). For variants with significant marginal genetic effects, SPAGxE<sub>CCT</sub> additionally calculates p value through Wald test and uses Cauchy combination (CCT) to combine p values from Wald test and the proposed genotype-adjusted test statistics.



# Introduction of Retrospective saddlepoint approximation (Retrospective SPA) approach in GWAS

Retrospective saddlepoint approximation (Retrospective SPA) is a method applied in GWAS. For a score statistic (S=G^TR), retrospective SPA strategy considers genotypes as random variables and approximate the distribution of score statistics conditional on phenotype and covariates.

Suggestions or comments on retrospective saddlepoint approximation methods are also welcome.

# SPAmix can control for population admixture

In the SPAmix paper (**A scalable, accurate, and universal analysis framework using individual-level allele frequency for large-scale genetic association studies in an admixed population**), we proposed SPAmix framework that is applicable to admixed populations and can incorporate local ancestry information into analyses. 

