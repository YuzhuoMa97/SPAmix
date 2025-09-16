# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/power/hetero_effect_size_scen1/
# sbatch -J QQplot --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript power_plot_v1.R $SLURM_ARRAY_TASK_ID'


library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(stringr)

sigLevel = 5e-8


MAF_ance1_List = c(0.01, 0.1)               # MAF in ancestry 1
MAF_ance2_List = c(0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3) # (causal ancestry)

# prev_List = c(0.2, 0.3)
prev_List = c(0.2) # 0.3 can not be analyzed by GRAB SPAmix
rep_List = c(1:100)

params = expand.grid(prev = prev_List,
                     MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     Gamma_ance2 = Gamma_ance2_Vec)


power = lapply(1 : nrow(params),function(i){
  print(i)
  MAF_ance1 = params[i, "MAF_ance1"]
  MAF_ance2 = params[i, "MAF_ance2"]
  Gamma_ance2 = params[i, "Gamma_ance2"]
  prev = params[i, "prev"]
  
  cat("\tMAF_ance1:\t",MAF_ance1,
      "\tMAF_ance2:\t",MAF_ance2,
      "\tGamma_ance2:\t",Gamma_ance2,
      "\tprevlence:\t",prev, "\n")
  
  files_SPAmix_index = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen1/prev_",prev,
                                         "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance_causal_",Gamma_ance2),pattern = ".index",full.names = T)
  file.remove(files_SPAmix_index)   ## remove index
  
  file_tractor = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen1/prev_",prev,
                                   "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance_causal_",Gamma_ance2),full.names = T) %>% str_subset("tsv")
  
  file_SPAmix = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/hetero_effect_size_scen1/prev_",prev,
                                  "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance_causal_",Gamma_ance2),full.names = T) %>% str_subset("SPAmix")
  
  results_tractor_temp = lapply(file_tractor, function(f){
    read.csv(f, sep="\t") 
  }) %>% do.call("rbind",.)%>%as_tibble()
  
  results_SPAmix_temp = lapply(file_SPAmix, function(f){
    fread(f) 
  }) %>% do.call("rbind",.)%>%as_tibble()
  
  results_temp = data.frame(prev = paste0("Prev = ", prev),
                            MAF_ance1 = paste0("MAF in ance2 = ", MAF_ance1),
                            MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
                            Gamma_ance2 = Gamma_ance2,
                            SPAmix_Pvalue = results_SPAmix_temp$Pvalue,
                            Tractor_ance1_Pvalue = results_tractor_temp$Gpval_anc0,
                            Tractor_ance2_Pvalue = results_tractor_temp$Gpval_anc1
  )
  
  power_SPAmix = nrow(results_temp %>% filter(SPAmix_Pvalue < sigLevel))/ nrow(results_temp)
  power_Tractor_ance1 = nrow(results_temp %>% filter(Tractor_ance1_Pvalue < sigLevel))/ nrow(results_temp)
  power_Tractor_ance2 = nrow(results_temp %>% filter(Tractor_ance2_Pvalue < sigLevel))/ nrow(results_temp)
  power_Tractor = nrow(results_temp %>% filter(Tractor_ance1_Pvalue < sigLevel | Tractor_ance2_Pvalue < sigLevel))/ nrow(results_temp)
  
  power_temp = data.frame(prev = paste0("Prev = ", prev),
                          MAF_ance1 = paste0("MAF in ance1 = ", MAF_ance1),
                          MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
                          Gamma_ance2 = Gamma_ance2,
                          SPAmix = power_SPAmix,
                          Tractor_ance1_null_ance = power_Tractor_ance1,
                          Tractor_ance2_causal_ance = power_Tractor_ance2,
                          Tractor = power_Tractor)
  
}) %>%do.call("rbind",.) %>% as_tibble() 

power = as.data.frame(power)

powerdata = reshape2::melt(power, id = c("prev", "MAF_ance1", "MAF_ance2", "Gamma_ance2"))%>%
  rename(Method = variable, power = value)

powerdata$Method <- factor(powerdata$Method,
                           levels = c("SPAmix",
                                      "Tractor_ance1_null_ance",
                                      "Tractor_ance2_causal_ance",
                                      "Tractor"))

powerdata = powerdata %>% arrange(Method)

# p1 = ggplot(powerdata, aes(x = Gamma_ance1, y = power, colour = Method, shape = Method, group = Method,linetype = Method)) +
#   geom_point(alpha = 0.7, size = 1.1) +
#   geom_line(linewidth = 0.7) + 
#   facet_grid(MAF_ance1 ~ MAF_ance2) +
#   labs(x="Genetic Effect", y="Empirical Power", title = paste0("Prevalence = ", prev_List[1])) +
#   theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

p1 = ggplot(powerdata, aes(x = Gamma_ance2, y = power, colour = Method, shape = Method, group = Method,linetype = Method)) +
  geom_point(alpha = 0.7, size = 1.25) +
  geom_line(linewidth = 0.7) + 
  facet_grid(MAF_ance1 ~ MAF_ance2) +
  # labs(x="Genetic Effect", y="Empirical Power", title = paste0("Prevalence = ", prev_List[1])) +
  labs(x="True Genetic Effect Size of Causal Ancestry", y="Empirical Powers") +
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

# ggsave(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/fig/power_prev_",prev_List[1],".jpeg"), p1)

ggsave(p1, filename = paste("power_prev",prev_List[1],".jpeg"), width = 8.5, height = 4.5,
       dpi = 400,
       path = "/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/fig/hetero_effect_size_scen1/")

