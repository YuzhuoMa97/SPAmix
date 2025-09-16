# cd /gdata01/user/yuzhuoma/SPA-G/UKB/all_population/code/
# sbatch --exclude=node05 -J step1 --mem=5000M -t 20-0:0 --array=1-12 -o log/%A_%a.log --wrap='Rscript Phenotype_SPAmix_SPAmixnoSPA_SPACox_SPACoxWB_allchr_step1_v2.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)  

library(GRAB)
library(data.table)
library(SPACox)
library(dplyr)
library(tidyr)
library(reshape2)

load("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/data/SampleIDs_all.RData") # 338,041 SampleIDs_all
load("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/data/SampleIDs_WB.RData")  # 281,296 SampleIDs_WB

PheCodeVec = list.files("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/")

# for (PheCode in PheCodeVec) {
PheCode = PheCodeVec[n.cpu]
PhenoFile = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/",PheCode),pattern = ".csv")
Phenotype = unlist(strsplit(PhenoFile, "_TTE.csv"))
cat("PheCode:\t", PheCode, "\n")
cat("PhenoFile:\t", PhenoFile, "\n")
cat("Phenotype:\t", Phenotype, "\n")
# Phenotype = "AD"

load("/gdata01/user/home/wenjianb/LA_PD/SPAGE_2022-06-23/version-1/before-step2-table.RData")
PhenoData = data.table::fread(paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/",PheCode, "/", PhenoFile)) # phenotype for all population

# PhenoData = data.table::fread("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/X290.11/AD_TTE.csv") # phenotype for all population
wjb = data.table::fread("/gdata02/master_data1/UK_Biobank/Data_Tools/dataMergeUmich_2022-05-09.csv")
PheData = merge(PhenoData, wjb)

##phenotypes in the population 
PheData_all = PheData %>% filter(Wenjian %in% as.numeric(SampleIDs_all)) %>%
  mutate(time = age.end) %>%
  select(Wenjian, birthYear, sex_genetic, 
         PC_ALL_1, PC_ALL_2, PC_ALL_3, PC_ALL_4, PC_ALL_5,
         PC_ALL_6, PC_ALL_7, PC_ALL_8, PC_ALL_9, PC_ALL_10,
         time, event, age.start, age.end, if_WB) %>%
  rename(IID = Wenjian) %>%
  arrange(IID)

PheData_WB = PheData_all %>% filter(IID %in% as.numeric(SampleIDs_WB)) %>% arrange(IID)


###fitting the null  model ---------------------------------------------------------------



objNull_SPAmix_All = GRAB.NullModel(Surv(age.start,age.end,event) ~ PC_ALL_1 + PC_ALL_2 + PC_ALL_3+ PC_ALL_4 + 
                                      PC_ALL_5 + PC_ALL_6 + PC_ALL_7 + PC_ALL_8 + PC_ALL_9 + PC_ALL_10 + sex_genetic + birthYear,
                                    data = PheData_all,
                                    subjData = IID,
                                    method = "SPAmix",
                                    traitType = "time-to-event",
                                    control = list(PC_columns = c('PC_ALL_1,PC_ALL_2,PC_ALL_3,PC_ALL_4,PC_ALL_5,PC_ALL_6,PC_ALL_7,PC_ALL_8,PC_ALL_9,PC_ALL_10')))

objNull_SPAmix_WB = GRAB.NullModel(Surv(age.start,age.end,event) ~ PC_ALL_1 + PC_ALL_2 + PC_ALL_3+ PC_ALL_4 + 
                                     PC_ALL_5 + PC_ALL_6 + PC_ALL_7 + PC_ALL_8 + PC_ALL_9 + PC_ALL_10 + sex_genetic + birthYear,
                                   data = PheData_WB,
                                   subjData = IID,
                                   method = "SPAmix",
                                   traitType = "time-to-event",
                                   control = list(PC_columns = c('PC_ALL_1,PC_ALL_2,PC_ALL_3,PC_ALL_4,PC_ALL_5,PC_ALL_6,PC_ALL_7,PC_ALL_8,PC_ALL_9,PC_ALL_10')))

save(objNull_SPAmix_All, file = paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/",PheCode, "/", "objNull_SPAmix_All_v2.RData"))
save(objNull_SPAmix_WB, file = paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/",PheCode, "/", "objNull_SPAmix_WB_v2.RData"))

# }




