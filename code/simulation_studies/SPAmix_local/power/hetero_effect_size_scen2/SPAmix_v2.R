# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/power/hetero_effect_size_scen2/
# sbatch --exclude=node05 -J SPAmix --mem=4000M -t 1-0:0 --array=1-1032 -o log/%A_%a.log --wrap='Rscript SPAmix_v2.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
library(data.table)
library(GRAB)

# Gets NULL model residuals for SPA-G
SPA_G_Get_Resid = function(traits="survival/binary",
                           formula=NULL,
                           data=NULL,
                           pIDs=NULL,
                           gIDs=NULL,
                           range=c(-100,100),
                           length.out = 10000,
                           ...)
  
{
  if(traits=="survival"){
    Call = match.call()
    ### Fit a Cox model
    obj.coxph = coxph(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.coxph, range)
    
    ### Get the covariate matrix to adjust for genotype
    resid = obj.coxph$residuals
    
    re = resid
  }
  else if(traits=="binary"){
    Call = match.call()
    ### Fit a logistic model
    obj.logistic = glm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.logistic, range)
    
    ### Get the covariate matrix to adjust for genotype
    mu = obj.logistic$fitted.values
    resid = obj.logistic$y - mu
    re = resid
  }
  return(re)
}

check_input = function(pIDs, gIDs, obj, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = obj$residuals
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}

check_input_Resid = function(pIDs, gIDs, R, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = R
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}



###################################################################################

MAF_ance1_List = c(0.01, 0.1)               # MAF in ancestry 1
MAF_ance2_List = c(0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
# Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3) # non-fixed ancestry 2

# Gamma_ance2_Vec = -c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3) # non-fixed ancestry 2
# Gamma_ance1_Vec = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 2.5) # Fixed ancestry 1

Gamma_ance1_Vec = c(0, 0.5, 1) # Fixed ancestry 1
Gamma_ance2_Vec = c(seq(0, 2, 0.05), 2.5, 3)   # non-fixed ancestry 2

prev_List = c(0.2)
rep_List = c(1:100)



params = expand.grid(prev = prev_List,
                     MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     Gamma_ance1 = Gamma_ance1_Vec,
                     Gamma_ance2 = Gamma_ance2_Vec)


MAF_ance1 = params[n.cpu, "MAF_ance1"] # MAF in ancestry 1
MAF_ance2 = params[n.cpu, "MAF_ance2"] # MAF in ancestry 2
prev = params[n.cpu, "prev"] #prevalence
rep = params[n.cpu, "rep"] # rep

MAF_ance1 = params[n.cpu, "MAF_ance1"] # MAF in ancestry 1
MAF_ance2 = params[n.cpu, "MAF_ance2"] # MAF in ancestry 2
prev = params[n.cpu, "prev"] #prevalence
gamma_ance1 = params[n.cpu, "Gamma_ance1"] # Fix genetic effect size of ancestry 1 
gamma_ance2 = params[n.cpu, "Gamma_ance2"] # genetic effect size of non-fixed ancestry 2 



dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output")
dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2")
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,"/gamma_ance2_",gamma_ance2))

#############################################################################################################
GenoFile = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/genotype/MAF_ance1_",
                  MAF_ance1, "_MAF_ance2_", MAF_ance2, "/Geno.mtx.bed")

load(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/genotype/MAF_ance1_",
            MAF_ance1, "_MAF_ance2_", MAF_ance2, "/Geno.mtx.RData"))


for (rep in rep_List) {
  print(rep)
  
  Pheno.mtx = data.table::fread(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/phenotype/hetero_effect_size_scen2/prev_",prev,
                                       "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,"/gamma_ance2_",gamma_ance2,
                                       "/PhenoData_rep",rep,".txt"),sep ="\t") %>% as.data.frame()
  
  
  ### fit null model
  resid  = SPA_G_Get_Resid(traits = "binary",
                           y~Cov1 + Cov2 + PC1 + PC2 + PC3 + PC4,family=binomial(link="logit"),
                           data=Pheno.mtx,
                           pIDs=Pheno.mtx$IID,
                           gIDs=rownames(Geno.mtx))
  
  objNull_SPAmix = GRAB.NullModel(resid ~ PC1 + PC2 + PC3 + PC4,
                                  data = Pheno.mtx, subjData = Pheno.mtx$IID, method = "SPAmix", traitType = "Residual", 
                                  control = list(PC_columns = "PC1,PC2,PC3,PC4"))
  
  #######
  
  
  
  output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                       "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,"/gamma_ance2_",gamma_ance2,
                       "/SPAmix_results_rep",rep,".txt")
  
  
  GRAB.Marker(objNull_SPAmix,
              GenoFile = GenoFile,
              OutputFile = output_file,
              control = list(AllMarkers = TRUE,
                             outputColumns = "zScore"))
  
  # GRAB.Marker(objNull_SPAmix,
  #             GenoFile = GenoFile,
  #             OutputFile = output_file,
  #             control = list(AllMarkers = TRUE,
  #                            AlleleOrder = "ref-first",
  #                            outputColumns = "zScore"))
  
  # GRAB.Marker(objNull_SPAmix,
  #             GenoFile = GenoFile,
  #             # GenoFileIndex = GenoFileIndex,
  #             OutputFile = output_file,
  #             control = list(AllMarkers = TRUE,
  #                            AlleleOrder = "ref-first",
  #                            min_maf_marker = 0.0001,
  #                            outputColumns = "zScore"))
}

