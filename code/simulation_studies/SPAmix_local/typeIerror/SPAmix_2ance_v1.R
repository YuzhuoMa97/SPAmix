# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/typeIerror/
# sbatch --exclude=node05 -J SPAmix --mem=4000M -t 1-0:0 --array=1-3200 -o log/%A_%a.log --wrap='Rscript SPAmix_2ance_v1.R $SLURM_ARRAY_TASK_ID'

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

source("/gdata01/user/yuzhuoma/SPA-G/tractor/code/Cauchy_combination_test.R")

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
prev_List = c(0.01, 0.05, 0.1, 0.2)
rep_List = c(1:100)



params = expand.grid(MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     prev = prev_List,
                     rep = rep_List)


MAF_ance1 = params[n.cpu, "MAF_ance1"] # MAF in ancestry 1
MAF_ance2 = params[n.cpu, "MAF_ance2"] # MAF in ancestry 2
prev = params[n.cpu, "prev"] #prevalence
rep = params[n.cpu, "rep"] # rep


dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                  MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev),recursive = T)



load(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
            MAF_ance1, "_MAF_ance2_", MAF_ance2, "/Geno.mtx.RData"))

#############################################################################################################


load(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
            MAF_ance1, "_MAF_ance2_", MAF_ance2, "/Geno.mtx.RData"))

load(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
            MAF_ance1, "_MAF_ance2_", MAF_ance2, "/Geno.mtx.ance1.RData"))

load(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
            MAF_ance1, "_MAF_ance2_", MAF_ance2, "/Geno.mtx.ance2.RData"))

haplo1.raw = data.table::fread(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
                                      MAF_ance1, "_MAF_ance2_", MAF_ance2, "/MAF_ance1_",MAF_ance1,"_MAF_ance2_",MAF_ance2,".anc1.hapcount.txt"))
haplo2.raw = data.table::fread(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
                                      MAF_ance1, "_MAF_ance2_", MAF_ance2, "/MAF_ance1_",MAF_ance1,"_MAF_ance2_",MAF_ance2,".anc2.hapcount.txt"))


haplo1 = haplo1.raw[,-1:-5]
haplo2 = haplo2.raw[,-1:-5]

MAF_ance1_est = colSums(Geno.mtx.ance1) / rowSums(haplo1)
MAF_ance2_est = colSums(Geno.mtx.ance2) / rowSums(haplo2)


# for (rep in rep_List) {
print(rep)

Pheno.mtx = data.table::fread(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/phenotype/prev_",prev,
                                     "/PhenoData_rep",rep,".txt"),sep ="\t") %>% as.data.frame()


### fit null model
resid  = SPA_G_Get_Resid(traits = "binary",
                         y~Cov1 + Cov2 + PC1 + PC2 + PC3 + PC4,family=binomial(link="logit"),
                         data=Pheno.mtx,
                         pIDs=Pheno.mtx$IID,
                         gIDs=rownames(Geno.mtx))

# objNull_SPAmix = GRAB.NullModel(resid ~ PC1 + PC2 + PC3 + PC4,
#                                 data = Pheno.mtx, subjData = Pheno.mtx$IID, method = "SPAmix", traitType = "Residual", 
#                                 control = list(PC_columns = "PC1,PC2,PC3,PC4"))

#######
pvalTot = c()

for(i in 1:ncol(Geno.mtx))
{
  print(i)
  # weight = dbeta(MAF_est[[i]], par1, par2)
  # wResid = weight * Data$Resid
  # print(mean(wResid))
  # wResid = wResid - mean(wResid)
  
  # Score.old = sum(resid * Geno.mtx[[i]])
  # Score.new = sum(wResid * Geno[[i]])
  Score.new1 = sum(resid * Geno.mtx.ance1[,i])
  Score.new2 = sum(resid * Geno.mtx.ance2[,i])
  
  # Var.old = sum(resid^2 * 2 * MAF_est[[i]] * (1-MAF_est[[i]]))
  
  mean.new1 = MAF_ance1_est[i] * sum(as.numeric(haplo1[i,])*resid)
  mean.new2 = MAF_ance2_est[i] * sum(as.numeric(haplo2[i,])*resid)
  
  Var.new1 = sum(resid^2 * as.numeric(haplo1[i,]) * MAF_ance1_est[i] * (1-MAF_ance1_est[i]))
  Var.new2 = sum(resid^2 * as.numeric(haplo2[i,]) * MAF_ance2_est[i] * (1-MAF_ance2_est[i]))
  
  # z.old = Score.old / sqrt(Var.old)
  z.new1 = (Score.new1 - mean.new1) / sqrt(Var.new1)
  z.new2 = (Score.new2 - mean.new2) / sqrt(Var.new2)
  
  # pval.old = 2*pnorm(abs(z.old), lower.tail=FALSE)
  pval.SPAmix.ance1 = 2*pnorm(abs(z.new1), lower.tail=FALSE)
  pval.SPAmix.ance2 = 2*pnorm(abs(z.new2), lower.tail=FALSE)
  
  # pval_CTT = CCT(c(pval.SPAmix.ance1, pval.SPAmix.ance2))
  
  pvalTot = rbind(pvalTot, c(colnames(Geno.mtx)[i], pval.SPAmix.ance1, pval.SPAmix.ance2))
  
  # pvalTot = rbind(pvalTot, c(colnames(Geno.mtx)[i], pval.SPAmix.ance1, pval.SPAmix.ance2, pval_CTT))
}

pvalTot = as.data.frame(pvalTot)
# colnames(pvalTot) = c("Marker", "pval.SPAmix.ance1", "pval.SPAmix.ance2", "pval_CTT")
colnames(pvalTot) = c("Marker", "pval.SPAmix.ance1", "pval.SPAmix.ance2")

output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                     MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev, "/SPAmix_2ance_results_rep",rep,".txt")

data.table::fwrite(pvalTot, file = output_file)
# }

