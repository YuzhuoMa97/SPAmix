# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/power/homo_effect_size/
# sbatch --partition=bi1 --exclude=node06 -J tractor --mem=4000M -t 1-0:0 --array=1-328 -o log/%A_%a.log --wrap='Rscript tractor_v1.R $SLURM_ARRAY_TASK_ID'

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

MAF_ance1_List = c(0.01, 0.1)               # MAF in ancestry 1
MAF_ance2_List = c(0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
# Gamma_ance1_Vec = Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
#                                       1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2, 2.5, 3)

Gamma_ance1_Vec = Gamma_ance2_Vec = seq(0, 2, 0.05)

prev_List = c(0.2)
rep_List = c(1:100)

params = expand.grid(prev = prev_List,
                     MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     Gamma_ance1 = Gamma_ance1_Vec)
                     # rep = rep_List)


MAF_ance1 = params[n.cpu, "MAF_ance1"] # MAF in ancestry 1
MAF_ance2 = params[n.cpu, "MAF_ance2"] # MAF in ancestry 2
prev = params[n.cpu, "prev"] #prevalence
gamma_ance1 = params[n.cpu, "Gamma_ance1"] # genetic effect size of ancestry 1
gamma_ance2 = gamma_ance1 # genetic effect size of ancestry 2
# rep = params[n.cpu, "rep"] # rep


dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/MAF_ance1_",
                  MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev),recursive = T)

dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output")
dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size")
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1))


hapdose_file =  paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/genotype/MAF_ance1_",
                       MAF_ance1, "_MAF_ance2_", MAF_ance2, "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2)

for (rep in rep_List) {
  print(rep)
  
  pheno_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/phenotype/homo_effect_size/prev_",prev,
                      "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,
                      "/PhenoData_rep",rep,".txt")
  
  output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                       "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,
                       "/tractor_results_rep",rep,".tsv")
  
  system(paste0("python /gdata01/user/yuzhuoma/Tractor/Tractor-master/RunTractor.py --hapdose ", 
                hapdose_file,
                " --phe ",pheno_file, " --method logistic --out ", output_file))
}






