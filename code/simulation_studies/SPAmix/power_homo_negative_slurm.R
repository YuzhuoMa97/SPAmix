# cd /gdata01/user/yuzhuoma/SPA-G/power/Rcode/5e-8/GRAB_SPAmix/
# sbatch --exclude=node05,node02 -J powerhomo --mem=4000M -t 1-0:0 --array=1-900 -o log/%A_%a.log --wrap='Rscript power_homo_negative_slurm.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

#### 3 Methods: SPACox, SPAmix, SPAmix-NoSPA, and add MAC and MAF.est.negative.number
library(data.table)
library(ggplot2)
library(dplyr)
library(survival)
library(SPACox)
library(GRAB)
# source("/gdata01/user/yuzhuoma/SPA-G/Rcode/SPAmix-functions-simulation-2023-02-23.R")
# source("/gdata01/user/yuzhuoma/SPA-G/Rcode/SPAmix-functions-cpp-2023-03-13.R")

# --------------------------------------------------------------------------------
N = 10000  # sample size
load("/gdata01/user/yuzhuoma/SPA-G/data/PCdata/top10PCs-1e5SNPs-MAFsumover0.1/top10PCs_MAFsumover0.1.RData")
load("/gdata01/user/yuzhuoma/SPA-G/data/a.population.structure.6.RData")

a = a6  # ancestry vector from EUR and EAS
PC = top10PCs # top 10 PCs of admixed populations calculated using 1e6 common variants

#### top 10 PCs
PC1 = PC[,1]
PC2 = PC[,2]
PC3 = PC[,3]
PC4 = PC[,4]
PC5 = PC[,5]
PC6 = PC[,6]
PC7 = PC[,7]
PC8 = PC[,8]
PC9 = PC[,9]
PC10 = PC[,10]

top4PCs = cbind(PC1, PC2, PC3, PC4)

# used in SPA-G.PC
top4PCs.intercept = cbind(1,top4PCs) # used in SPA-G.PC
PC.invPCPC = top4PCs.intercept %*% solve(t(top4PCs.intercept)%*%top4PCs.intercept) # used in SPA-G.PC
tPC = t(top4PCs.intercept) # used in SPA-G.PC

#### Genotype matrix
#### min.MAF is the minimal MAF of some SNP in EUR and EAS

load("/gdata01/user/yuzhuoma/SPA-G/power/data/genotype-MAFsumover0.002/Geno.mtx1.RData")  # min.MAF < 0.01
load("/gdata01/user/yuzhuoma/SPA-G/power/data/genotype-MAFsumover0.002/Geno.mtx2.RData")  # 0.01 =< min.MAF < 0.05
load("/gdata01/user/yuzhuoma/SPA-G/power/data/genotype-MAFsumover0.002/Geno.mtx3.RData")  # min.MAF >= 0.05

Geno.mtx1.file = "/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/genotype-MAFsumover0.002/Geno.mtx1.bed" # min.MAF < 0.01
Geno.mtx2.file = "/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/genotype-MAFsumover0.002/Geno.mtx2.bed" # 0.01 =< min.MAF < 0.05
Geno.mtx3.file = "/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/genotype-MAFsumover0.002/Geno.mtx3.bed" # min.MAF >= 0.05

#-----------------------------------------------------------------------
nSNP = dim(Geno.mtx1)[2]/5  # SNP numbers = 1000
Pheno.number = 1000          # number of phenotype
Pheno.merge.number = 10        # merge Pheno.merge.num phenoptypes in one n.cpu
Pheno.scenraio.number = Pheno.number/Pheno.merge.number

#### Scenarios
par.ls = c()
for(ER in c(1, 2, 3)){
  for (min.MAF in c(1, 2, 3)) {
    for (Pheno.scenraio.ID in 1:Pheno.scenraio.number) {
      par.ls = rbind(par.ls, c(ER, min.MAF, Pheno.scenraio.ID))
    }
  }
}

colnames(par.ls) = c("ER", "min.MAF", "Pheno.scenraio.ID")
# n.cpu = 900
ER = par.ls[n.cpu,"ER"]
min.MAF = par.ls[n.cpu,"min.MAF"]

Pheno.scenraio.ID = par.ls[n.cpu,"Pheno.scenraio.ID"]

##### PRINT
ER
min.MAF
Pheno.scenraio.ID

#### Event rate of ancestry populations (EUR and EAS)
ER1Vec = c(0.01, 0.05, 0.2)   # orginal event rate for subpopulation 1 (EUR)
ER2Vec = c(0.01, 0.05, 0.2)  # orginal event rate for subpopulation 2 (EAS)

if(ER == 1){
  ER1 = ER1Vec[1]
  ER2 = ER2Vec[1]
  ER_level = "L_ER"
}
if(ER == 2){
  ER1 = ER1Vec[2]
  ER2 = ER2Vec[2]
  ER_level = "M_ER"
}
if(ER == 3){
  ER1 = ER1Vec[3]
  ER2 = ER2Vec[3]
  ER_level = "H_ER"
}


#### minimal of MAF in EUR and EAS
if(min.MAF == 1){
  Minimal.AF = "Minimal.AF < 0.01"
  minMAF_level = "L_minMAF"
  Geno.mtx = Geno.mtx1
  GenoFile = Geno.mtx1.file
}

if(min.MAF == 2){
  Minimal.AF = "0.01 <= Minimal.AF < 0.05"
  minMAF_level = "M_minMAF"
  Geno.mtx = Geno.mtx2
  GenoFile = Geno.mtx2.file
}

if(min.MAF == 3){
  Minimal.AF = "Minimal.AF >= 0.05"
  minMAF_level = "H_minMAF"
  Geno.mtx = Geno.mtx3
  GenoFile = Geno.mtx3.file
}

####  list of scenarios of MAF difference in EUR and EAS, and each Geno.mtx include the following five scenarios, that is, 
####  the first thousand columns of Geno.mtx are genotype of SNPs whose "MAF.EUR - MAF.EAS <= -0.05"
####  the sceond thousand columns of Geno.mtx are genotype of SNPs whose "-0.05 < MAF.EUR - MAF.EAS <= -0.01"
####  the third thousand columns of Geno.mtx are genotype of SNPs whose "-0.01 < MAF.EUR - MAF.EAS < 0.01"
####  the fourth thousand columns of Geno.mtx are genotype of SNPs whose "0.01 <= MAF.EUR - MAF.EAS < 0.05"
####  the fifth thousand columns of Geno.mtx are genotype of SNPs whose "MAF.EUR - MAF.EAS >= 0.05"

table.Diff = data.table::data.table(Diff = c(rep("DiffMAF << 0",nSNP),
                                             rep("DiffMAF < 0",nSNP),
                                             rep("DiffMAF ~ 0",nSNP),
                                             rep("DiffMAF > 0",nSNP), 
                                             rep("DiffMAF >> 0",nSNP)))

# dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/power/data/power-v19-20230418-log10MAF/Pheno.scenraio.ID-",Pheno.scenraio.ID),recursive = T)
# file.output = paste0("/gdata01/user/yuzhuoma/SPA-G/power/data/power-v19-20230418-log10MAF/Pheno.scenraio.ID-",Pheno.scenraio.ID,"/min.MAF-",min.MAF,
#                      "-ER.EUR-",ER1,"-ER.EAS-",ER2,".csv")

dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power_homo_negative/GRAB_SPAmix/Pheno.scenraio.ID-",Pheno.scenraio.ID),recursive = T)
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power_homo_negative/SPACox/Pheno.scenraio.ID-",Pheno.scenraio.ID),recursive = T)

# file.output = paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power-GRAB_SPAmix-2log10MAF/Pheno.scenraio.ID-",Pheno.scenraio.ID,"/min.MAF-",min.MAF,
#                      "-ER.EUR-",ER1,"-ER.EAS-",ER2,".csv")

#### Calculate p values ################################################################################################
# pvalMat = c() # prepare matrix of p-values


for (Pheno.ID in c((Pheno.scenraio.ID-1)*Pheno.merge.number + 1:Pheno.merge.number)) {
  
  # phenotype data include true ancestry proportions of EAS and two covariates of X1(standard normal distribution) and X2(Bernoulli distribution)
  load(paste0("/gdata01/user/yuzhuoma/SPA-G/power/data/phenotype-homo-MAFsumover0.002-2log10MAF-neggamma/event.rate-",ER,
              "/Phenotype-", Pheno.ID,".Rdata"))
  
  
  Phen.mat = data.frame(ID = paste0("Subj-", 1:10000),
                        event = Phen.mtx$event,
                        time = Phen.mtx$time,
                        Cov1 = PC1,
                        Cov2 = PC2,
                        Cov3 = PC3,
                        Cov4 = PC4,
                        Cov5 = PC5,
                        Cov6 = PC6,
                        Cov7 = PC7,
                        Cov8 = PC8,
                        Cov9 = PC9,
                        Cov10 = PC10,
                        Cov11 = Phen.mtx$X2,
                        Cov12 = Phen.mtx$X3)
  
  #### calculate martingale residuals
  objNull_SPAmix = GRAB.NullModel(Surv(time,event)~Cov1+Cov2+Cov3+Cov4+Cov11+Cov12,
                                  data=Phen.mat,
                                  subjData=Phen.mat$ID,
                                  method = "SPAmix",
                                  traitType = "time-to-event",
                                  control = list(PC_columns = c('Cov1,Cov2,Cov3,Cov4')))
  
  obj.null.SPACox = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2+Cov3+Cov4+Cov11+Cov12,
                                      data=Phen.mat,
                                      pIDs=Phen.mtx$ID,
                                      gIDs=paste0("IID-",1:N))
  
  
  #### calculate p-values from different methods
  
  # OutputFile = paste0("/gdata01/user/yuzhuoma/SPA-G/data/Type-one-error-rate_GRAB_SPAmix_1e6/Pheno.scenraio.ID-",Pheno.scenraio.ID,".txt")
  dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power_homo_negative/GRAB_SPAmix/Pheno.scenraio.ID-",Pheno.scenraio.ID,"/Pheno.ID-",Pheno.ID),recursive = T)
  dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power_homo_negative/SPACox/Pheno.scenraio.ID-",Pheno.scenraio.ID,"/Pheno.ID-",Pheno.ID),recursive = T)
  
  file.output.SPAmix = paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power_homo_negative/GRAB_SPAmix/Pheno.scenraio.ID-",
                              Pheno.scenraio.ID,"/Pheno.ID-",Pheno.ID,"/min.MAF-",min.MAF,"-ER.EUR-",ER1,"-ER.EAS-",ER2,".txt")
  
  GRAB.Marker(objNull_SPAmix,
              GenoFile = GenoFile,
              OutputFile = file.output.SPAmix,
              control = list(AllMarkers = TRUE,
                             min_maf_marker = 0.0001,
                             outputColumns = "zScore",
                             min_mac_marker = 1))
  
  
  
  
  
  # SPACox
  file.output.SPACox = paste0("/gdata01/user/yuzhuoma/SPA-G/power/5e-8/data/power_homo_negative/SPACox/Pheno.scenraio.ID-",
                              Pheno.scenraio.ID,"/Pheno.ID-",Pheno.ID,"/min.MAF-",min.MAF,"-ER.EUR-",ER1,"-ER.EAS-",ER2,".csv")
  
  res.SPACox = SPACox(obj.null = obj.null.SPACox, Geno.mtx = Geno.mtx, min.maf = 0.00001)
  
  
  pvalMat.SPACox =  cbind(table.Diff, Method = "SPACox", res.SPACox[,"p.value.spa"]) %>% as.data.frame()
  colnames(pvalMat.SPACox) = c("DiffMAF", "Method", "pvalue")
  rownames(pvalMat.SPACox) = c(1:nrow(pvalMat.SPACox))
  pvalMat.SPACox$pvalue = as.numeric(pvalMat.SPACox$pvalue)
  
  pvalMat.SPACox$DiffMAF <- factor(pvalMat.SPACox$DiffMAF,levels=c("DiffMAF << 0",
                                                                   "DiffMAF < 0",
                                                                   "DiffMAF ~ 0",
                                                                   "DiffMAF > 0",
                                                                   "DiffMAF >> 0"))
  
  pvalMat.SPACox = cbind("Minimal.AF" = Minimal.AF,
                         "minMAF_level" = minMAF_level,
                         "ER" = paste0("ER in EUR = ", ER1, sep=", ", "ER in EAS = ", ER2),
                         "ER_level" = ER_level,
                         pvalMat.SPACox)
  
  pvalMat.SPACox = as.data.frame(pvalMat.SPACox)
  
  data.table::fwrite(pvalMat.SPACox, file.output.SPACox)
  
}

