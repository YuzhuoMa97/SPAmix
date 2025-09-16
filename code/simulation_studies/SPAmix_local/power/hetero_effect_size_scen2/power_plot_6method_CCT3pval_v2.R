# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/power/hetero_effect_size_scen2/
# sbatch -J QQplot --mem=4000M -t 1-0:0 --array=1-3 -o log/%A_%a.log --wrap='Rscript power_plot_6method_CCT3pval_v2.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(stringr)

source("/gdata01/user/yuzhuoma/SPA-G/tractor/code/Cauchy_combination_test.R")

sigLevel = 5e-8


MAF_ance1_List = c(0.01, 0.1)               # MAF in ancestry 1
MAF_ance2_List = c(0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
# Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3) # non-fixed ancestry 2
# Gamma_ance2_Vec_pos = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3) # non-fixed ancestry 2
# Gamma_ance2_Vec_neg = -c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3) # non-fixed ancestry 2
# Gamma_ance2_Vec = c(Gamma_ance2_Vec_neg, Gamma_ance2_Vec_pos)

Gamma_ance2_Vec = c(seq(0, 2, 0.05), 2.5, 3)   # non-fixed ancestry 2
Gamma_ance1_Vec = c(0, 0.5, 1) # Fixed ancestry 1


# Gamma_ance1_Vec = c(0, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 2.5) # Fixed ancestry 1
# Gamma_ance1_Vec = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 2.5) # Fixed ancestry 1
prev_List = c(0.2) # 0.3 can not be analyzed by GRAB SPAmix
rep_List = c(1:100)

params = expand.grid(prev = prev_List,
                     MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     Gamma_ance2 = Gamma_ance2_Vec)

Gamma_ance1 = Gamma_ance1_Vec[n.cpu] # Fix genetic effect size of ancestry 1 



power = lapply(1 : nrow(params),function(i){
  print(i)
  MAF_ance1 = params[i, "MAF_ance1"]
  MAF_ance2 = params[i, "MAF_ance2"]
  Gamma_ance2 = params[i, "Gamma_ance2"]
  prev = params[i, "prev"]
  
  cat("\tMAF_ance1:\t",MAF_ance1,
      "\tMAF_ance2:\t",MAF_ance2,
      "\tGamma_ance2:\t",Gamma_ance2,
      "\tprevalence:\t",prev, "\n")
  
  files_SPAmix_index = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                                         "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1,"/gamma_ance2_",Gamma_ance2),pattern = ".index",full.names = T)
  file.remove(files_SPAmix_index)   ## remove index
  
  file_tractor = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                                   "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1,"/gamma_ance2_",Gamma_ance2),full.names = T) %>% str_subset("tsv")
  
  # # file of results of GRAB SPAmix and SPAmix_2ance (return p-values by ancestry)
  # file_SPAmix_0 = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
  #                                   "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1,"/gamma_ance2_",Gamma_ance2),full.names = T) %>% str_subset("SPAmix")
  # 
  # # file of results of SPAmix_2ance (return p-values by ancestry)
  # file_SPAmix_2ance = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
  #                                       "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1,"/gamma_ance2_",Gamma_ance2),full.names = T) %>% str_subset("2ance")
  
  # file of results of GRAB SPAmix, SPAmix_2ance, and SPAmix_CCT (combine p-values pg GRAB SPAmix, SPAmix_ance1 and SPAmix_ance2) (return p-values by ancestry)
  file_SPAmix_CCT = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen2/prev_",prev,
                                      "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1,"/gamma_ance2_",Gamma_ance2),full.names = T) %>% str_subset("CCT")
  
  
  results_tractor_temp = lapply(file_tractor, function(f){
    read.csv(f, sep="\t") 
  }) %>% do.call("rbind",.)%>%as_tibble()
  
  results_SPAmix_temp = lapply(file_SPAmix_CCT, function(f){
    fread(f) 
  }) %>% do.call("rbind",.)%>%as_tibble()
  

  results_temp = data.frame(prev = paste0("Prev = ", prev),
                            MAF_ance1 = paste0("MAF in ance2 = ", MAF_ance1),
                            MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
                            Gamma_ance2 = Gamma_ance2,
                            SPAmix_Pvalue = results_SPAmix_temp$Pvalue.GRAB.SPAmix,
                            SPAmix_ance1_Pvalue = results_SPAmix_temp$Pvalue.SPAmix.ance1,
                            SPAmix_ance2_Pvalue = results_SPAmix_temp$Pvalue.SPAmix.ance2,
                            SPAmix_CCT_Pvalue = results_SPAmix_temp$Pvalue_SPAmix_CCT,
                            Tractor_ance1_Pvalue = results_tractor_temp$Gpval_anc0,
                            Tractor_ance2_Pvalue = results_tractor_temp$Gpval_anc1
  )
  
  power_SPAmix = nrow(results_temp %>% filter(SPAmix_Pvalue < sigLevel))/ nrow(results_temp)
  power_SPAmix_ance1 = nrow(results_temp %>% filter(SPAmix_ance1_Pvalue < sigLevel))/ nrow(results_temp)
  power_SPAmix_ance2 = nrow(results_temp %>% filter(SPAmix_ance2_Pvalue < sigLevel))/ nrow(results_temp)
  power_SPAmix_CCT = nrow(results_temp %>% filter(SPAmix_CCT_Pvalue < sigLevel))/ nrow(results_temp)
  
  
  power_Tractor_ance1 = nrow(results_temp %>% filter(Tractor_ance1_Pvalue < sigLevel))/ nrow(results_temp)
  power_Tractor_ance2 = nrow(results_temp %>% filter(Tractor_ance2_Pvalue < sigLevel))/ nrow(results_temp)

  # power_Tractor = nrow(results_temp %>% filter(Tractor_ance1_Pvalue < sigLevel | Tractor_ance2_Pvalue < sigLevel))/ nrow(results_temp)
  
  power_temp = data.frame(prev = paste0("Prev = ", prev),
                          MAF_ance1 = paste0("MAF in ance1 = ", MAF_ance1),
                          MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
                          Gamma_ance2 = Gamma_ance2,
                          SPAmix = power_SPAmix,
                          SPAmix_ance1_fixed = power_SPAmix_ance1,
                          SPAmix_ance2_nonfixed = power_SPAmix_ance2,
                          SPAmix_CCT = power_SPAmix_CCT,
                          Tractor_ance1_fixed = power_Tractor_ance1,
                          Tractor_ance2_nonfixed = power_Tractor_ance2)
  
}) %>%do.call("rbind",.) %>% as_tibble() 

power = as.data.frame(power)

powerdata = reshape2::melt(power, id = c("prev", "MAF_ance1", "MAF_ance2", "Gamma_ance2"))%>%
  rename(Method = variable, power = value)

powerdata$Method <- factor(powerdata$Method,
                           levels = c("SPAmix",
                                      "SPAmix_ance1_fixed",
                                      "SPAmix_ance2_nonfixed",
                                      "SPAmix_CCT",
                                      "Tractor_ance1_fixed",
                                      "Tractor_ance2_nonfixed"))

powerdata = powerdata %>% arrange(Method)

# p1 = ggplot(powerdata, aes(x = Gamma_ance1, y = power, colour = Method, shape = Method, group = Method,linetype = Method)) +
#   geom_point(alpha = 0.7, size = 1.1) +
#   geom_line(linewidth = 0.7) + 
#   facet_grid(MAF_ance1 ~ MAF_ance2) +
#   labs(x="Genetic Effect", y="Empirical Power", title = paste0("Prevalence = ", prev_List[1])) +
#   theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

p1 = ggplot(powerdata, aes(x = Gamma_ance2, y = power, colour = Method, shape = Method, group = Method,linetype = Method)) +
  geom_point(alpha = 0.7, size = 1.25) +
  geom_vline(xintercept = Gamma_ance1, linetype = "dashed", color = "grey", linewidth = 0.5, alpha = 0.8)+
  geom_line(linewidth = 0.7) + 
  facet_grid(MAF_ance1 ~ MAF_ance2) +
  # labs(x="Genetic Effect", y="Empirical Power", title = paste0("Prevalence = ", prev_List[1])) +
  labs(x="True Genetic Effect Size of Ancestry 2", y="Empirical Powers") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom") +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.title = element_blank(),
    # legend.position="none",
    # panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank()
    # panel.grid.minor = element_blank()
    # legend.title=element_blank()
  ) 
# theme(strip.text.x = element_text(size = 12),
#       strip.text.y = element_text(size = 12))+ 
# theme(axis.title.x = element_text(size = 14),
#       axis.title.y = element_text(size = 14)) +
# theme(legend.text = element_text(size=13),
#       legend.title = element_text(size=15))


dir.create("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/fig/hetero_effect_size_scen2/6Methods_comb/")

ggsave(p1, filename = paste("power_prev",prev_List[1],"_gamma_ancestry1_",Gamma_ance1,".jpeg"), width = 12, height = 4.5,
       dpi = 400,
       path = "/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/fig/hetero_effect_size_scen2/6Methods_comb/")

data.table::fwrite(powerdata, file = paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/fig/hetero_effect_size_scen2/6Methods_comb/power_prev",prev_List[1],"gamma_ancestry1_",Gamma_ance1,".csv"))