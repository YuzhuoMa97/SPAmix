# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/typeIerror/
# sbatch --exclude=node05 -J SPAmix --mem=4000M -t 1-0:0 --array=1-3200 -o log/%A_%a.log --wrap='Rscript SPAmix_CCT_v1.R $SLURM_ARRAY_TASK_ID'

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



print(rep)


output_file_GRAB_SPAmix = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                                 MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev, "/SPAmix_results_rep",rep,".txt")


output_file_SPAmix_2ance =  paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                                   MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_",prev, "/SPAmix_SPA_2ance_results_rep",rep,".txt")



data_GRAB_SPAmix = data.table::fread(output_file_GRAB_SPAmix)
data_SPAmix_2ance = data.table::fread(output_file_SPAmix_2ance)

colnames(data_GRAB_SPAmix) = paste0(colnames(data_GRAB_SPAmix),".GRAB.SPAmix")
data_GRAB_SPAmix = data_GRAB_SPAmix %>% rename(Marker = Marker.GRAB.SPAmix)

data_temp = merge(data_GRAB_SPAmix, data_SPAmix_2ance)

results_SPAmix_CCT_temp = data_temp %>% mutate(Pvalue_SPAmix_CCT = NA)


#### add CTT p-values for tractor
# pval_SPAmix_CCT = c()
for (k in 1:nrow(data_temp)) {
  print(k)
  pval_SPAmix_CCT_temp = CCT(c(results_SPAmix_CCT_temp[k,]$Pvalue.GRAB.SPAmix,
                               results_SPAmix_CCT_temp[k,]$Pvalue.SPAmix.ance1,
                               results_SPAmix_CCT_temp[k,]$Pvalue.SPAmix.ance2))
  
  results_SPAmix_CCT_temp$Pvalue_SPAmix_CCT[k] = pval_SPAmix_CCT_temp
}

# results_SPAmix_CCT_temp = data_temp %>% mutate(pval_SPAmix_CCT = pval_SPAmix_CCT)


output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/output/MAF_ance1_",
                     MAF_ance1, "_MAF_ance2_", MAF_ance2, "_prev_", prev, "/SPAmix_CCT_3pval_results_rep", rep,".txt")

data.table::fwrite(results_SPAmix_CCT_temp, file = output_file)

