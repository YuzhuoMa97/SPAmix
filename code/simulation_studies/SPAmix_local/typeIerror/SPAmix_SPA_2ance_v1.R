# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/typeIerror/
# sbatch --exclude=node05 -J SPAmix --mem=4000M -t 1-0:0 --array=1-3200 -o log/%A_%a.log --wrap='Rscript SPAmix_SPA_2ance_v1.R $SLURM_ARRAY_TASK_ID'

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
source("/gdata01/user/yuzhuoma/SPA-G/tractor/code/SPAmix-functions-local-ance-2023-08-28.R")

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

haplo.mtx.ance1 = t(haplo1)
haplo.mtx.ance2 = t(haplo2)

colnames(haplo.mtx.ance1) = colnames(Geno.mtx.ance1)
colnames(haplo.mtx.ance2) = colnames(Geno.mtx.ance2)

# MAF_ance1_est = colSums(Geno.mtx.ance1) / rowSums(haplo1)
# MAF_ance2_est = colSums(Geno.mtx.ance2) / rowSums(haplo2)


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
#######
pvalTot_ance1 = SPAmix_localance(Geno.mtx = Geno.mtx.ance1,
                                 R = resid,
                                 haplo.mtx = haplo.mtx.ance1,
                                 # obj.null,
                                 Cutoff = 2,
                                 impute.method = "fixed",
                                 missing.cutoff = 0.15,
                                 min.maf = 0.00001,          
                                 G.model = "Add")

colnames(pvalTot_ance1) = c("Marker", "MAF.ance1","missing.rate.ance1","Pvalue.SPAmix.ance1","Pvalue.norm.ance1",
                            "Stat.ance1","Mean.ance1","Var.ance1","z.ance1", "MAC.ance1")


pvalTot_ance2 = SPAmix_localance(Geno.mtx = Geno.mtx.ance2,
                                 R = resid,
                                 haplo.mtx = haplo.mtx.ance2,
                                 # obj.null,
                                 Cutoff = 2,
                                 impute.method = "fixed",
                                 missing.cutoff = 0.15,
                                 min.maf = 0.00001,          
                                 G.model = "Add")

colnames(pvalTot_ance2) = c("Marker", "MAF.ance2","missing.rate.ance2","Pvalue.SPAmix.ance2","Pvalue.norm.ance2",
                            "Stat.ance2","Mean.ance2","Var.ance2","z.ance2", "MAC.ance2")

pvalTot_ance1_ance2 = merge(pvalTot_ance1, pvalTot_ance2)


output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                     MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev, "/SPAmix_SPA_2ance_results_rep",rep,".txt")

data.table::fwrite(pvalTot_ance1_ance2, file = output_file)


