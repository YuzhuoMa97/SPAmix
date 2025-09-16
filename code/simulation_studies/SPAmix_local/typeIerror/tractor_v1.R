# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/typeIerror/
# sbatch --exclude=node06 -J tractor --mem=4000M -t 1-0:0 --array=1-1600 -o log/%A_%a.log --wrap='Rscript tractor_v1.R $SLURM_ARRAY_TASK_ID'

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

# # MAF_ance1_List = c(0.001, 0.01, 0.1)               # MAF in ancestry 1
# # MAF_ance1_List = c(0.001, 0.01)               # MAF in ancestry 1
# MAF_ance1_List = c(0.1)               # MAF in ancestry 1
# 
# MAF_ance2_List = c(0.001, 0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
# 
# prev_List = c(0.001, 0.01, 0.05, 0.3)
# rep_List = c(1:100)


MAF_ance1_List = c(0.01, 0.1)               # MAF in ancestry 1

MAF_ance2_List = c(0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
prev_List = c(0.01, 0.05)

# prev_List = c(0.1, 0.2)
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

# setwd(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
#              MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev))


hapdose_file =  paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/genotype/MAF_ance1_",
                       MAF_ance1, "_MAF_ance2_", MAF_ance2, "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2)

pheno_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/phenotype/prev_",prev,
                    "/PhenoData_rep",rep,".txt")

output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                     MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev, "/tractor_results_rep",rep,".tsv")

system(paste0("python /gdata01/user/yuzhuoma/Tractor/Tractor-master/RunTractor.py --hapdose ", 
              hapdose_file,
              " --phe ",pheno_file, " --method logistic --out ", output_file))


# test = read.csv("Y:/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_0.001_MAF_ance2_0.001_prev_0.001/results_rep1.tsv",
#                 sep = "\t")


