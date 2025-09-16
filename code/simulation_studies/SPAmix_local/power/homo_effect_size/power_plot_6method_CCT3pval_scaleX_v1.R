# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/power/homo_effect_size/
# sbatch -J QQplot --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript power_plot_6method_CCT3pval_scaleX_v1.R $SLURM_ARRAY_TASK_ID'


library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(stringr)

sigLevel = 5e-8


MAF_ance1_List = c(0.01, 0.1)               # MAF in ancestry 1
MAF_ance2_List = c(0.01, 0.05, 0.1, 0.3)    # MAF in ancestry 2
# Gamma_ance1_Vec = Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3)
# Gamma_ance1_Vec = Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5)
# Gamma_ance1_Vec = Gamma_ance2_Vec = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
#                                       1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2)


Gamma_ance1_Vec = Gamma_ance2_Vec = seq(0, 2, 0.05)
prev_List = c(0.2) # 0.3 can not be analyzed by GRAB SPAmix
rep_List = c(1:100)

params = expand.grid(prev = prev_List,
                     MAF_ance1 = MAF_ance1_List,
                     MAF_ance2 = MAF_ance2_List,
                     Gamma_ance1 = Gamma_ance1_Vec)


power = lapply(1 : nrow(params),function(i){
  print(i)
  MAF_ance1 = params[i, "MAF_ance1"]
  MAF_ance2 = params[i, "MAF_ance2"]
  Gamma_ance1 = params[i, "Gamma_ance1"]
  prev = params[i, "prev"]
  
  cat("\tMAF_ance1:\t",MAF_ance1,
      "\tMAF_ance2:\t",MAF_ance2,
      "\tGamma_ance1:\t",Gamma_ance1,
      "\tprevlence:\t",prev, "\n")
  
  files_SPAmix_index = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                                         "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1),pattern = ".index",full.names = T)
  file.remove(files_SPAmix_index)   ## remove index
  
  file_tractor = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                                   "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1),full.names = T) %>% str_subset("tsv")
 
  # file of results of GRAB SPAmix, SPAmix_2ance, and SPAmix_CCT (combine p-values pg GRAB SPAmix, SPAmix_ance1 and SPAmix_ance2) (return p-values by ancestry)
  file_SPAmix_CCT = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/output/homo_effect_size/prev_",prev,
                                      "/MAF_ance1_", MAF_ance1, "_MAF_ance2_", MAF_ance2,"/gamma_ance1_",Gamma_ance1),full.names = T) %>% str_subset("CCT")
  
  
  results_tractor_temp = lapply(file_tractor, function(f){
    read.csv(f, sep="\t") 
  }) %>% do.call("rbind",.)%>%as_tibble()
  
  results_SPAmix_temp = lapply(file_SPAmix_CCT, function(f){
    fread(f) 
  }) %>% do.call("rbind",.)%>%as_tibble()
  

  
  results_temp = data.frame(prev = paste0("Prev = ", prev),
                            MAF_ance1 = paste0("MAF in ance2 = ", MAF_ance1),
                            MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
                            Gamma_ance1 = Gamma_ance1,
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
  # 
  # power_temp = data.frame(prev = paste0("Prev = ", prev),
  #                         MAF_ance1 = paste0("MAF in ance1 = ", MAF_ance1),
  #                         MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
  #                         Gamma_ance1 = Gamma_ance1,
  #                         SPAmix = power_SPAmix,
  #                         SPAmix_ance1_fixed = power_SPAmix_ance1,
  #                         SPAmix_ance2_nonfixed = power_SPAmix_ance2,
  #                         SPAmix_CCT = power_SPAmix_CCT,
  #                         Tractor_ance1_fixed = power_Tractor_ance1,
  #                         Tractor_ance2_nonfixed = power_Tractor_ance2)
  
  
  power_temp = data.frame(prev = paste0("Prev = ", prev),
                          MAF_ance1 = paste0("MAF in ance1 = ", MAF_ance1),
                          MAF_ance2 = paste0("MAF in ance2 = ", MAF_ance2),
                          Gamma_ance1 = Gamma_ance1,
                          SPAmix = power_SPAmix,
                          SPAmix_ance1 = power_SPAmix_ance1,
                          SPAmix_ance2 = power_SPAmix_ance2,
                          SPAmix_CCT = power_SPAmix_CCT,
                          Tractor_ance1 = power_Tractor_ance1,
                          Tractor_ance2 = power_Tractor_ance2)
  
  
}) %>%do.call("rbind",.) %>% as_tibble() 

power = as.data.frame(power)

powerdata = reshape2::melt(power, id = c("prev", "MAF_ance1", "MAF_ance2", "Gamma_ance1"))%>%
  rename(Method = variable, power = value)

powerdata$Method <- factor(powerdata$Method,
                           levels = c("SPAmix",
                                      "SPAmix_ance1",
                                      "SPAmix_ance2",
                                      "SPAmix_CCT",
                                      "Tractor_ance1",
                                      "Tractor_ance2"))

powerdata = powerdata %>% arrange(Method)

#################################################################################################################

power_MAFance1_0.01 = powerdata %>% filter(MAF_ance1 == "MAF in ance1 = 0.01")
power_MAFance1_0.01_MAFance2_0.01 = power_MAFance1_0.01 %>% filter(MAF_ance2 == "MAF in ance2 = 0.01")
power_MAFance1_0.01_MAFance2_0.05 = power_MAFance1_0.01 %>% filter(MAF_ance2 == "MAF in ance2 = 0.05")
power_MAFance1_0.01_MAFance2_0.1 = power_MAFance1_0.01 %>% filter(MAF_ance2 == "MAF in ance2 = 0.1")
power_MAFance1_0.01_MAFance2_0.3 = power_MAFance1_0.01 %>% filter(MAF_ance2 == "MAF in ance2 = 0.3")

power_MAFance1_0.01_MAFance2_0.01_plot = power_MAFance1_0.01_MAFance2_0.01  %>%
  filter(Gamma_ance1 %in% seq(0, 2, 0.2))

power_MAFance1_0.01_MAFance2_0.05_plot = power_MAFance1_0.01_MAFance2_0.05  %>%
  filter(Gamma_ance1 %in% seq(0, 1, 0.1))

power_MAFance1_0.01_MAFance2_0.1_plot = power_MAFance1_0.01_MAFance2_0.1  %>%
  filter(Gamma_ance1 %in% seq(0, 0.8, 0.05))

power_MAFance1_0.01_MAFance2_0.3_plot = power_MAFance1_0.01_MAFance2_0.3 %>%
  filter(Gamma_ance1 %in% seq(0, 0.5, 0.05))

power_MAFance1_0.01_plot = rbind.data.frame(power_MAFance1_0.01_MAFance2_0.01_plot,
                                            power_MAFance1_0.01_MAFance2_0.05_plot,
                                            power_MAFance1_0.01_MAFance2_0.1_plot,
                                            power_MAFance1_0.01_MAFance2_0.3_plot)


fig1 = ggplot(power_MAFance1_0.01_plot, aes(x = Gamma_ance1, y = power, colour = Method, shape = Method, group = Method,linetype = Method)) +
  geom_point(alpha = 0.7, size = 1.25) +
  geom_line(linewidth = 0.7) + 
  facet_grid(MAF_ance1 ~ MAF_ance2, scales = "free_x") +
  labs(x="", y="Empirical Powers") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank()
  ) 


power_MAFance1_0.1 = powerdata %>% filter(MAF_ance1 == "MAF in ance1 = 0.1")
power_MAFance1_0.1_MAFance2_0.01 = power_MAFance1_0.1 %>% filter(MAF_ance2 == "MAF in ance2 = 0.01")
power_MAFance1_0.1_MAFance2_0.05 = power_MAFance1_0.1 %>% filter(MAF_ance2 == "MAF in ance2 = 0.05")
power_MAFance1_0.1_MAFance2_0.1 = power_MAFance1_0.1 %>% filter(MAF_ance2 == "MAF in ance2 = 0.1")
power_MAFance1_0.1_MAFance2_0.3 = power_MAFance1_0.1 %>% filter(MAF_ance2 == "MAF in ance2 = 0.3")

power_MAFance1_0.1_MAFance2_0.01_plot = power_MAFance1_0.1_MAFance2_0.01  %>%
  filter(Gamma_ance1 %in% seq(0, 0.7, 0.05))

power_MAFance1_0.1_MAFance2_0.05_plot = power_MAFance1_0.1_MAFance2_0.05  %>%
  filter(Gamma_ance1 %in% seq(0, 0.7, 0.05))

power_MAFance1_0.1_MAFance2_0.1_plot = power_MAFance1_0.1_MAFance2_0.1  %>%
  filter(Gamma_ance1 %in% seq(0, 0.7, 0.05))

power_MAFance1_0.1_MAFance2_0.3_plot = power_MAFance1_0.1_MAFance2_0.3 %>%
  filter(Gamma_ance1 %in% seq(0, 0.7, 0.05))

power_MAFance1_0.1_plot = rbind.data.frame(power_MAFance1_0.1_MAFance2_0.01_plot,
                                           power_MAFance1_0.1_MAFance2_0.05_plot,
                                           power_MAFance1_0.1_MAFance2_0.1_plot,
                                           power_MAFance1_0.1_MAFance2_0.3_plot)



fig2 = ggplot(power_MAFance1_0.1_plot, aes(x = Gamma_ance1, y = power, colour = Method, shape = Method, group = Method,linetype = Method)) +
  geom_point(alpha = 0.7, size = 1.25) +
  geom_line(linewidth = 0.7) + 
  facet_grid(MAF_ance1 ~ MAF_ance2, scales = "free_x") +
  labs(x="True Genetic Effect Size", y="Empirical Powers") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom") +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank()
  ) 



# library(cowplot)
# power_fig_part1 = plot_grid(fig1, fig2, fig3, fig4, fig5, ncol = 1)
# power_fig_part2 = fig6
library(gridExtra)
power_fig = gridExtra::grid.arrange(fig1, fig2, heights=c(1.87, 2.5))

ggsave(power_fig, filename = paste("power_prev",prev_List[1],"_scaleX.jpeg"), width = 10, height = 5.5,
       dpi = 400,
       path = "/gdata01/user/yuzhuoma/SPA-G/tractor/data/power_v1/fig/homo_effect_size")

