
library(dplyr)
library(tidyr)

load("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/topPCs/top10PCs.RData")
top10PCs = as.data.frame(top10PCs)

load("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/global_ancestry/global_ancestry.RData")
gloval_ance = data.frame(Ancestry1 = 1-alpha,
                         Ancestry2 = alpha)

library(ggplot2)
QQplot1 = ggplot(top10PCs, aes(Comp.1, Comp.2)) +
  geom_point(alpha = 0.3) + labs(x="PC1",y="PC2") + theme_bw() 

QQplot2 = ggplot(top10PCs, aes(Comp.2, Comp.3)) +
  geom_point(alpha = 0.3) + labs(x="PC2",y="PC3") + theme_bw() 

data_ancestryVec_PCs = cbind.data.frame(gloval_ance, top10PCs)

QQplot3 = ggplot(data_ancestryVec_PCs, aes(Ancestry1, Comp.1)) +
  geom_point(alpha = 0.3) + labs(x="Ancestry Proportion of Ancestry1",y="PC1") + theme_bw() 


gloval_ance_v2 = gloval_ance %>% arrange(Ancestry1) %>% mutate(index = c(1:10000))
gloval_ance_v3 = rbind(cbind(gloval_ance_v2[,"Ancestry1"],gloval_ance_v2[,"index"],"Ancestry 1"),
                       cbind(gloval_ance_v2[,"Ancestry2"],gloval_ance_v2[,"index"],"Ancestry 2")) %>%
  as.data.frame() %>% rename(Global_ance = V1, index = V2, Ancestry = V3) %>%
  mutate(Global_ance = as.numeric(Global_ance), index = as.numeric(index))



Fig4 = ggplot(data = gloval_ance_v3, aes(x=index,y=Global_ance, fill = Ancestry)) +
  geom_bar(stat = "identity") + theme(plot.margin = margin(0)) + labs(x="Individuals",y="Global Ancestry Proportions") +
  theme_bw() + theme(legend.position = "bottom") +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # strip.background = element_blank(),
    # legend.text.align = 0  # Left-align the legend text
  ) 

library(cowplot)
QQplot = plot_grid(QQplot1, QQplot2, QQplot3, Fig4,
                   ncol = 2, labels = c("a.", "b.", "c.", "d."))




ggsave(QQplot, filename ="ance_PCfig.jpeg",
       width = 9, height = 8,
       dpi = 400,
       path = "/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/fig/")











