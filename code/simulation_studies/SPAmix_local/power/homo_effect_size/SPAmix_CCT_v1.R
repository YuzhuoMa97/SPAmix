# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/power/homo_effect_size/
# sbatch --partition=bi1 --exclude=node05 -J SPAmix --mem=4000M -t 1-0:0 --array=1-328 -o log/%A_%a.log --wrap='Rscript SPAmix_CCT_v1.R $SLURM_ARRAY_TASK_ID'

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
# Gamma_ance1_Vec = Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
#                                       1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2, 2.5, 3)

Gamma_ance1_Vec = Gamma_ance2_Vec = seq(0, 2, 0.05)


# prev_List = c(0.2, 0.3)
prev_List = c(0.2) # 0.3 can not be analyzed by GRAB SPAmix
rep_List = c(1:100)



params = expand.grid(prev = prev_List,
                     MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     Gamma_ance1 = Gamma_ance1_Vec)


MAF_ance1 = params[n.cpu, "MAF_ance1"] # MAF in ancestry 1
MAF_ance2 = params[n.cpu, "MAF_ance2"] # MAF in ancestry 2
prev = params[n.cpu, "prev"] #prevalence
rep = params[n.cpu, "rep"] # rep

MAF_ance1 = params[n.cpu, "MAF_ance1"] # MAF in ancestry 1
MAF_ance2 = params[n.cpu, "MAF_ance2"] # MAF in ancestry 2
prev = params[n.cpu, "prev"] #prevalence
gamma_ance1 = params[n.cpu, "Gamma_ance1"] # genetic effect size of ancestry 1
gamma_ance2 = gamma_ance1 # genetic effect size of ancestry 2



dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output")
dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size")
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2))
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1))

#############################################################################################################
for (rep in rep_List) {
  print(rep)
  
  output_file_GRAB_SPAmix = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                                   "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,
                                   "/SPAmix_results_rep",rep,".txt")
  
  # output_file_SPAmix_2ance = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
  #                                   "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,"/gamma_ance2_",gamma_ance2,
  #                                   "/SPAmix_2ance_results_rep",rep,".txt")
  
  output_file_SPAmix_2ance = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                                    "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,
                                    "/SPAmix_SPA_2ance_results_rep",rep,".txt")
  
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
  
  
  
  output_file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                       "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",gamma_ance1,
                       "/SPAmix_CCT_3pval_results_rep",rep,".txt")
  
  data.table::fwrite(results_SPAmix_CCT_temp, file = output_file)
}



