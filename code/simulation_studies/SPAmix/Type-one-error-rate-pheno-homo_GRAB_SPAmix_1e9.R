# cd /gdata01/user/yuzhuoma/SPA-G/Rcode/
# sbatch --exclude=node05,node06,node04 -J 9homo --mem=5000M -t 30-0:0 --array=1-90 -o log/%A_%a.log --wrap='Rscript Type-one-error-rate-pheno-homo_GRAB_SPAmix_1e9.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

#### 2 Methods: SPAmix and SPAmix-NoSPA (use GRAB package)
library(data.table)
library(GRAB)
library(survival)

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

#### Genotype matrix (AF.EUR + AF.EAS > 0.002)
#### min.MAF is the minimal MAF of some SNP in EUR and EAS
# load("/gdata01/user/yuzhuoma/SPA-G/data/genotype/Geno.mtx1.RData")  # min.MAF < 0.01
# load("/gdata01/user/yuzhuoma/SPA-G/data/genotype/Geno.mtx2.RData")  # 0.01 =< min.MAF < 0.05
# load("/gdata01/user/yuzhuoma/SPA-G/data/genotype/Geno.mtx3.RData")  # min.MAF >= 0.05

Geno.mtx1 = "/gdata01/user/yuzhuoma/SPA-G/data/genotype/Geno.mtx1.bed" # min.MAF < 0.01
Geno.mtx2 = "/gdata01/user/yuzhuoma/SPA-G/data/genotype/Geno.mtx2.bed" # 0.01 =< min.MAF < 0.05
Geno.mtx3 = "/gdata01/user/yuzhuoma/SPA-G/data/genotype/Geno.mtx3.bed" # min.MAF >= 0.05

#-----------------------------------------------------------------------
# nSNP = dim(Geno.mtx1)[2]/5  # SNP numbers = 1000
nSNP = 1000                      # SNP numbers = 1000
Pheno.number = 1000000              # number of phenotype
Pheno.merge.number = 100000        # merge Pheno.merge.num phenoptypes in one n.cpu
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
# n.cpu = 1
ER = par.ls[n.cpu,"ER"]
min.MAF = par.ls[n.cpu,"min.MAF"]

Pheno.scenraio.ID = par.ls[n.cpu,"Pheno.scenraio.ID"]

##### PRINT
ER
min.MAF
Pheno.scenraio.ID

#### Event rate of ancestry populations (EUR and EAS)
ER1Vec = c(0.01, 0.05, 0.2)    # orginal event rate for subpopulation 1 (EUR)
ER2Vec = c(0.01, 0.05, 0.2)    # orginal event rate for subpopulation 2 (EAS)

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
  GenoFile = Geno.mtx1
}

if(min.MAF == 2){
  Minimal.AF = "0.01 <= Minimal.AF < 0.05"
  minMAF_level = "M_minMAF"
  GenoFile = Geno.mtx2
}

if(min.MAF == 3){
  Minimal.AF = "Minimal.AF >= 0.05"
  minMAF_level = "H_minMAF"
  GenoFile = Geno.mtx3
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

dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/data/Type-one-error-rate-pheno-homo_GRAB_SPAmix_1e9/Pheno.scenraio.ID-",Pheno.scenraio.ID),recursive = T)
# file.output = paste0("/gdata01/user/yuzhuoma/SPA-G/data/Type-one-error-rate_GRAB_SPAmix_1e6/Pheno.scenraio.ID-",Pheno.scenraio.ID,"/min.MAF-",min.MAF,
#                      "-ER.EUR-",ER1,"-ER.EAS-",ER2,".csv")

#### Calculate p values ################################################################################################

for (Pheno.ID in c((Pheno.scenraio.ID-1)*Pheno.merge.number + 1:Pheno.merge.number)) {
  
  # phenotype data include two covariates of X1(standard normal distribution) and X2(Bernoulli distribution)
  load(paste0("/gdata01/user/yuzhuoma/SPA-G/data/phenotype-homo-1e6/event.rate-",ER,
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
  
  #### calculate p-values from SPAmix
  
  # OutputFile = paste0("/gdata01/user/yuzhuoma/SPA-G/data/Type-one-error-rate_GRAB_SPAmix_1e6/Pheno.scenraio.ID-",Pheno.scenraio.ID,".txt")
  dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/data/Type-one-error-rate-pheno-homo_GRAB_SPAmix_1e9/Pheno.scenraio.ID-",Pheno.scenraio.ID,"/Pheno.ID-",Pheno.ID),recursive = T)
  
  file.output = paste0("/gdata01/user/yuzhuoma/SPA-G/data/Type-one-error-rate-pheno-homo_GRAB_SPAmix_1e9/Pheno.scenraio.ID-",
                       Pheno.scenraio.ID,"/Pheno.ID-",Pheno.ID,"/min.MAF-",min.MAF,"-ER.EUR-",ER1,"-ER.EAS-",ER2,".txt")
  
  GRAB.Marker(objNull_SPAmix,
              GenoFile = GenoFile,
              OutputFile = file.output,
              control = list(AllMarkers = TRUE,
                             min_maf_marker = 0.0001,
                             outputColumns = "zScore",
                             min_mac_marker = 1))
  
}


