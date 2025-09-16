# cd /gdata01/user/yuzhuoma/SPA-G/UKB/all_population/code/
# sbatch --exclude=node06,node04,node05 -J Eczmix --mem=5000M -t 3-0:0 --array=1-2336 -o log/%A_%a.log --wrap='Rscript Phenotype_SPAmix_SPAmixnoSPA_allchr_step2.R $SLURM_ARRAY_TASK_ID'
# sbatch -J Parkinson --mem=7000M -t 3-0:0 --array=1-2336 -o log/%A_%a.log --wrap='Rscript Phenotype_SPAmix_SPAmixnoSPA_allchr_step2.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)  

library(GRAB)
library(data.table)

# updated on 2023-03-23: some variants sum(g.tilde)==0, thus MAF.est=0, and result in inflated type I error rates

# read in an object of "table"
load("/gdata01/user/home/wenjianb/LA_PD/SPAGE_2022-06-23/version-1/before-step2-table.RData")

# read in objNull from step 1
# PheCode = "X290.11" 
# PheCode = "X250.2" 
# PheCode = "X401.1" 
# PheCode = "X332" 
# PheCode = "X282.5" 
# PheCode = "X335" 
# PheCode = "X365" 
# PheCode = "X499" 
# PheCode = "X495" 
# PheCode = "X743.11" 
# PheCode = "X295.1" 
PheCode = "X931" 


PhenoFile = list.files(paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/",PheCode),pattern = ".csv")
Phenotype = unlist(strsplit(PhenoFile, "_TTE.csv"))
cat("PheCode:\t", PheCode, "\n")
cat("PhenoFile:\t", PhenoFile, "\n")
cat("Phenotype:\t", Phenotype, "\n")
# Phenotype = "AD" 
# load("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/X290.11/objNull_SPAmix_All.RData")
load(paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/working/",PheCode, "/", "objNull_SPAmix_All.RData"))

# Seperate SNPs into several parts
# n.cpu=1
files = table$file[table$group==n.cpu]
group = n.cpu
dir.create(paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/Results_All/", Phenotype, "/Results_SPAmix_SPAmixnoSPA"),recursive = T)
# file.output = paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/Results/", Phenotype, "_Results_SPAmix_SPAmixnoSPA/results_",n.cpu,".csv")

for(IDsToIncludeFile in files) {
# IDsToIncludeFile = files[3]
  IDsToIncludeFile = stringr::str_replace(IDsToIncludeFile, "/data1", "/gdata02/master_data1")
  cat("IDsToIncludeFile:\t",IDsToIncludeFile,"\n")
  chr = unlist(strsplit(IDsToIncludeFile, "chr_|_chunk"))[2]
  chunk = unlist(strsplit(IDsToIncludeFile, "chr_|_chunk_|.txt"))[3]
  chr_chunk = paste0("chr_",chr,"_chunk_",chunk)
  
  ##load in population data--------------------------------------------------------------
  ##Genotype data
  GenoFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c",chr,"_b0_v3.bgen")
  GenoFileIndex = c(paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr",chr,"_v3.bgen.bgi"),
                    
                    "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample");

  # OutputFile = paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/Results_All/", Phenotype, "/Results_SPAmix_SPAmixnoSPA/results_",n.cpu,".txt")
  OutputFile = paste0("/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/Results_All/", Phenotype, "/Results_SPAmix_SPAmixnoSPA/results_",chr_chunk,".txt")
  
  GRAB.Marker(objNull_SPAmix_All,
              GenoFile = GenoFile,
              GenoFileIndex = GenoFileIndex,
              OutputFile = OutputFile,
              control = list(IDsToIncludeFile = IDsToIncludeFile,
                             AlleleOrder = "ref-first",
                             min_maf_marker = 0.0001,
                             outputColumns = "zScore"))
  
}

# SampleIDs_all = intersect(rownames(GenoMat),as.character(PheData$Wenjian)) # 338,041
# save(SampleIDs_all, file = "/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/data/SampleIDs_all.RData")
# 
# SampleIDs_WB = intersect(rownames(GenoMat),as.character(PheData_WB$Wenjian)) # 281,296
# save(SampleIDs_WB, file = "/gdata01/user/yuzhuoma/SPA-G/UKB/all_population/data/SampleIDs_WB.RData")
