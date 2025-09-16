# ========================================================================
# SPAmix Analysis Pipeline for ALLofUS Database using GRAB Package
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Fit a Null Model for SPAmix Analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPAmix requires fitting a null model to account for covariates and population structure
# before performing marker-level association tests

# Mathematical Formulation:
# λ(t|X,G) = λ_0(t)exp(βX + γG)
# Where:
#   λ_0(t): Baseline hazard function
#   X: Fixed effects (age, genetic sex, PCs)

# -------------------------------
# 1.1 Data Preparation
# -------------------------------
# Load phenotype data (Replace with ALLofUS data path)
PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData = data.table::fread(PhenoFile, header = TRUE)

# Generate simulated principal components
# REAL DATA ANALYSIS: Replace with ALLofUS-provided genetic PCs
N = nrow(PhenoData)
PhenoData = PhenoData %>% mutate(
  PC1 = rnorm(N),  # Principal component 1
  PC2 = rnorm(N),  # Principal component 2
  PC3 = rnorm(N),  # Principal component 3
  PC4 = rnorm(N)   # Principal component 4
)

# ---------------------------------
# 1.2 Null Model Specification
# ---------------------------------
# Option 1: Direct survival analysis specification
# Best for: Initial analyses without pre-computed residuals
# Users can directly specify a time-to-event trait to analyze
obj.SPAmix = GRAB.NullModel(
  formula = Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2 + PC3 + PC4,  # a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (i.e. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis. Other values (e.g. -9, -999) will be treated as ordinary numeric values in analysis.
  data = PhenoData,                                                            # a data.frame, list or environment (or object coercible by as.data.frame to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
  subjData = IID,                                                              # Subject ID column (critical for data alignment): a character vector of subject IDs. Its order should be the same as the subject order in the formula and data (before subset process).
  method = "SPAmix", 
  traitType = "time-to-event",                                                 # SPAmix: Support traitType = "time-to-event" or "Residual". Please check GRAB.SPAmix for more details.
  control = list(                                                              # a list of parameters for controlling the model fitting process. For more details, please check Details section in GRAB.NullModel.
    PC_columns = "PC1,PC2,PC3,PC4"                                             # Must match PC names in data
  )
)

# Option 2: Residual-based approach (alternative specification)
# Useful when pre-computed residuals exist from external models
# Using model residuals performs exactly the same as the above. Note that confounding factors are still required in the right of the formula.

library(survival)
obj.coxph = coxph(
  Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2 + PC3 + PC4,
  data = PhenoData
)

# SPAmix null model using martingale residuals
obj.SPAmix = GRAB.NullModel(
  formula = obj.coxph$residuals ~ AGE + GENDER + PC1 + PC2 + PC3 + PC4,  # a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (i.e. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis. 
  data = PhenoData,                                                      # a data.frame, list or environment (or object coercible by as.data.frame to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
  subjData = IID,                                                        # Subject ID column (critical for data alignment): a character vector of subject IDs. Its order should be the same as the subject order in the formula and data (before subset process).
  method = "SPAmix",                                                     # a character
  traitType = "Residual",                                                # SPAmix: Support traitType = "time-to-event" or "Residual". Check GRAB.SPAmix for more details.
  control = list(                                                        # a list of parameters for controlling the model fitting process. For more details, please check Details section in GRAB.NullModel.
    PC_columns = "PC1,PC2,PC3,PC4"                                       # PC adjustment even with residuals (for calculating individual-specific allele frequencies)
  )
)

# Critical Validation Check:
# Please make sure the order of Subject ID in subjData is the same as the subject order in the formula, data and martingale residuals (obj.coxph$residuals).

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Genome-Wide Association Analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This section performs genetic association testing 

# -------------------------------
# File Path Configuration
# -------------------------------
# Please replace these with ALLofUS data paths
GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputDir = system.file("results", package = "GRAB")

# --------------------------------------------------
# 2A: Full Genome-Wide Analysis
# --------------------------------------------------
# Full genome analysis output file
OutputFile_AllMarkers_TRUE = paste0(OutputDir, "/Results_SPAmix_FullGWAS.txt")

GRAB.Marker(
  objNull = obj.SPAmix,                          # the output object of function GRAB.NullModel.
  GenoFile = GenoFile,                           # a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check GRAB.ReadGeno for more details.
  OutputFile = OutputFile_AllMarkers_TRUE,       # a character of output file to save the analysis results.
  control = list(                                # a list of parameters for controlling function GRAB.Marker, more details can be seen in Details section.
    nMarkersEachChunk = 10000,                   # number of markers (default=10000) in one chunk to output.
    AlleleOrder = "ref-first",                   # a character, "ref-first" or "alt-first", to determine whether the REF/major allele should appear first or second. Default is "alt-first" for PLINK and "ref-first" for BGEN. If the ALT allele frequencies of most markers are > 0.5, you should consider resetting this option. NOTE, if you use plink2 to convert PLINK file to BGEN file, then 'ref-first' modifier is to reset the order.
    AllMarkers = TRUE,                           # a logical value (default: FALSE) to indicate if all markers are extracted. It might take too much memory to put genotype of all markers in R. This parameter is to remind users.
    ImputeMethod = "mean",                       # a character, "none" (default), "bestguess", or "mean". By default, missing genotype is NA. Suppose alternative allele frequency is p, then missing genotype is imputed as 2p (ImputeMethod = "mean") or round(2p) (ImputeMethod = "bestguess").
    MissingRateCutoff = 0.15,                    # a numeric value (default=0.15). Markers with missing rate > this value will be excluded from analysis.
    min_maf_marker = 0.0001,                     # a numeric value (default=0.001). Markers with MAF < this value will be excluded from analysis.
    outputColumns = "zScore"                     # Output standardized test statistics
    # IDsToIncludeFile = my_IDsToIncludeFile,    # a file of marker IDs to include, one column (no header). Check system.file("extdata", "IDsToInclude.txt", package = "GRAB") for an example.
    # IDsToExcludeFile = IDsToExcludeFile,       # a file of marker IDs to exclude, one column (no header).
    # RangesToIncludeFile = RangesToIncludeFile, # a file of ranges to include, three columns (no headers): chromosome, start position, end position. Check system.file("extdata", "RangesToInclude.txt", package = "GRAB") for an example.
    # RangesToExcludeFile = RangesToExcludeFile, # a file of ranges to exclude, three columns (no headers): chromosome, start position, end position.
  )
)

# Quick results preview
data.table::fread(OutputFile_AllMarkers_TRUE)

# --------------------------------------------------
# 2B: Targeted Regional Analysis
# --------------------------------------------------
# Create inclusion list for targeted analysis 

# Get PLINK marker metadata
plink_file <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
geno_list <- GRAB.ReadGeno(
  plink_file,
  control = list(AllMarkers = TRUE)
)

# Generate inclusion file (first 5000 SNPs for demonstration)
first_5000_snps <- head(geno_list$markerInfo$ID, 5000)

# Write inclusion file (use tempfile() for production)
output_IDsToIncludeFile <- file.path(
  system.file(package = "GRAB"), 
  "extdata", 
  "Custom_SNP_List.txt"
)
write.table(
  first_5000_snps, 
  file = output_IDsToIncludeFile,
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE
)

# Perform targeted analysis
OutputFile_AllMarkers_FALSE = paste0(OutputDir, "/Results_SPAmix_Targeted.txt")

GRAB.Marker(
  objNull = obj.SPAmix,                          # the output object of function GRAB.NullModel.
  GenoFile = GenoFile,                           # a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check GRAB.ReadGeno for more details.
  OutputFile = OutputFile_AllMarkers_FALSE,      # a character of output file to save the analysis results.
  control = list(                                # a list of parameters for controlling function GRAB.Marker, more details can be seen in Details section.
    nMarkersEachChunk = 10000,                   # number of markers (default=10000) in one chunk to output.
    AlleleOrder = "ref-first",                   # a character, "ref-first" or "alt-first", to determine whether the REF/major allele should appear first or second. Default is "alt-first" for PLINK and "ref-first" for BGEN. If the ALT allele frequencies of most markers are > 0.5, you should consider resetting this option. NOTE, if you use plink2 to convert PLINK file to BGEN file, then 'ref-first' modifier is to reset the order.
    AllMarkers = TRUE,                           # a logical value (default: FALSE) to indicate if all markers are extracted. It might take too much memory to put genotype of all markers in R. This parameter is to remind users.
    ImputeMethod = "mean",                       # a character, "none" (default), "bestguess", or "mean". By default, missing genotype is NA. Suppose alternative allele frequency is p, then missing genotype is imputed as 2p (ImputeMethod = "mean") or round(2p) (ImputeMethod = "bestguess").
    MissingRateCutoff = 0.15,                    # a numeric value (default=0.15). Markers with missing rate > this value will be excluded from analysis.
    min_maf_marker = 0.0001,                     # a numeric value (default=0.001). Markers with MAF < this value will be excluded from analysis.
    outputColumns = "zScore",                    # Output standardized test statistics
    IDsToIncludeFile = output_IDsToIncludeFile   # a file of marker IDs to include, one column (no header). Check system.file("extdata", "IDsToInclude.txt", package = "GRAB") for an example.
    # IDsToExcludeFile = IDsToExcludeFile,       # a file of marker IDs to exclude, one column (no header).
    # RangesToIncludeFile = RangesToIncludeFile, # a file of ranges to include, three columns (no headers): chromosome, start position, end position. Check system.file("extdata", "RangesToInclude.txt", package = "GRAB") for an example.
    # RangesToExcludeFile                        # a file of ranges to exclude, three columns (no headers): chromosome, start position, end position.
  )
)

# Quick results preview
data.table::fread(OutputFile_AllMarkers_FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Post-Analysis Processing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Demonstration cleanup - omit for real analysis
file.remove(OutputFile_AllMarkers_TRUE)
file.remove(OutputFile_AllMarkers_FALSE)



