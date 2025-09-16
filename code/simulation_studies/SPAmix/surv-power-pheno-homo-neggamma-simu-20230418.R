# cd /gdata01/user/yuzhuoma/SPA-G/power/Rcode/
# sbatch --exclude=node05 -J jobID --mem=5000M -t 1-0:0 --array=1-3 -o log/%A_%a.log --wrap='Rscript surv-power-pheno-homo-neggamma-simu-20230418.R $SLURM_ARRAY_TASK_ID'
# sbatch -J phenosimu --mem=5000M -t 1-0:0 --array=1-3 -o log/%A_%a.log --wrap='Rscript surv-power-neggamma-pheno-homo-simu.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

N = 10000
N1 = 10000 # subpopulation 1 (EUR)
N2 = 10000 # subpopulation 2 (EAS)

load("/gdata01/user/yuzhuoma/a.population.structure.6.RData")
a = a6

# load genotype data
load("/gdata01/user/yuzhuoma/SPA-G/power/data/genotype-MAFsumover0.002/Geno.mtx.RData")
load("/gdata01/user/yuzhuoma/SPA-G/power/data/genotype-MAFsumover0.002/Geno.mtx.EUR.RData")
load("/gdata01/user/yuzhuoma/SPA-G/power/data/genotype-MAFsumover0.002/Geno.mtx.EAS.RData")

# Geno.mtx = as.matrix(Geno.mtx[,1])
# Geno.mtx.EUR = as.matrix(Geno.mtx.EUR[,1])
# Geno.mtx.EAS = as.matrix(Geno.mtx.EAS[,1])

Geno.mtx.adj = Geno.mtx
for (i in 1:N) {
  Geno.mtx.adj[i,] = Geno.mtx[i,] - colMeans(Geno.mtx)
}


# Geno.mtx.EUR.adj = Geno.mtx.EUR
# meanMAF = apply(Geno.mtx,2,mean)
# Geno.mtx.EUR.adj = apply(Geno.mtx.EUR.adj, 1,function(x){x-colMeans(Geno.mtx)})%>%t()

Geno.mtx.EUR.adj = Geno.mtx.EUR
for (i in 1:N) {
  Geno.mtx.EUR.adj[i,] = Geno.mtx.EUR[i,] - colMeans(Geno.mtx)
}


Geno.mtx.EAS.adj = Geno.mtx.EAS
for (i in 1:N) {
  Geno.mtx.EAS.adj[i,] = Geno.mtx.EAS[i,] - colMeans(Geno.mtx)
}

MAFVec = ifelse(colMeans(Geno.mtx)<0.5, colMeans(Geno.mtx), 1-colMeans(Geno.mtx))
gammaVec = 8*log10(MAFVec)

Pheno.number = 1000          # 1000 phenotypes
ER1Vec = c(0.01, 0.05, 0.2)    # orginal event rate for subpopulation 1 (EUR)
ER2Vec = c(0.01, 0.05, 0.2)    # orginal event rate for subpopulation 2 (EAS)

#####################################################################################################
library(survival)
#### lower function to estimate beta0 given an event rate. Will be used in data.simu.surv().
### simulate time-to-event phenotype for admixed population (e.g. EUR and EAS)

f1 = function(N1,
              event.rate1,
              lamda,
              gammaVec,
              g1.mtx)               # Genotype matrix of population 1
{
  scale0 = lamda
  shape0 = 2
  # set.seed(seed1)
  X21 = rnorm(N1)
  X31 = rbinom(N1, 1, 0.5)
  
  cens = rweibull(N1, shape=1, scale=0.15)
  eps <- runif(N1, 0, 1)
  time = (-log(eps)*exp(-X21*0.5-X31*0.5-g1.mtx%*%gammaVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  re = mean(event) - event.rate1
  return(re)
}

f2 = function(N1,
              N2,
              event.rate1,
              event.rate2,
              b,
              gammaVec,
              g1.mtx,       # Genotype matrix of population 1
              g2.mtx)       # Genotype matrix of population 2
  # seed1,
  # seed2)
{
  scale0 = uniroot(f1, c(-100000,100000), N1 = N1, event.rate1 = event.rate1, gammaVec = gammaVec, g1.mtx = g1.mtx)
  scale0 = scale0$root
  shape0 = 2
  # set.seed(seed2)
  # X1 = c(rep(0, N1), rep(1, N2))
  # X11 = X1[c(1 : N1)]
  X12 = rep(1, N2)
  X22 = rnorm(N2)
  X32 = rbinom(N2, 1, 0.5)
  
  cens = rweibull(N2, shape=1, scale=0.15)
  eps <- runif(N2, 0, 1)
  time = (-log(eps)*exp(-X12*b-X22*0.5-X32*0.5-g2.mtx%*%gammaVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  re = mean(event) - event.rate2
  return(re)
}

####  simulate time-to-event phenotype for population admixed with two subpopulations
####  event rates are different for different subjects
data.simu.surv.admix = function(N,
                                N1,
                                N2,
                                event.rate1,
                                event.rate2,
                                gammaVec,
                                g.mtx,               # genotype for N admixed population
                                g1.mtx,              # genotype for N EUR population
                                g2.mtx,              # genotype for N EAS population
                                # seed1,
                                # seed2,
                                a)
{
  scale0 = uniroot(f1, c(-100000,100000), N1 = N1, event.rate1 = event.rate1, gammaVec = gammaVec, g1.mtx = g1.mtx)
  scale0 = scale0$root
  shape0 = 2
  b = uniroot(f2, c(-100000,100000), N1 = N1, N2 = N2, event.rate1 = event.rate1, 
              event.rate2 = event.rate2, gammaVec = gammaVec, g1.mtx = g1.mtx, g2.mtx = g2.mtx)
  b = b$root
  
  # X1 = c(rep(0, N1), rep(1, N2))
  X2 = rnorm(N)
  X3 = rbinom(N, 1, 0.5)
  
  # set.seed(1)
  cens = rweibull(N, shape=1, scale=0.15) 
  eps <- runif(N, 0, 1)
  time = (-log(eps)*exp(-a[,2]*b-X2*0.5-X3*0.5-g.mtx%*%gammaVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  out = data.frame(surv.time, event, cens, a[,2], X2, X3)
  colnames(out)= c("surv.time","event", "cens", "X1.new", "X2", "X3")
  return(out)
}

###############################################################################
# gammaVec = runif(150, 0.5, 0.5)
# # gammaVec = c(runif(50, 2, 2), runif(50, 1, 1), runif(50, 0.5, 0.5)) 
gammaVec

ER = n.cpu

# for (ER in c(1, 2, 3)) {
if(ER == 1){
  ER1 = ER1Vec[1]
  ER2 = ER2Vec[1]
}
if(ER == 2){
  ER1 = ER1Vec[2]
  ER2 = ER2Vec[2]
}
if(ER == 3){
  ER1 = ER1Vec[3]
  ER2 = ER2Vec[3]
}

for (Pheno.ID in c(1:Pheno.number)) {
  cat(paste0("ER =", ER, sep=", ", "Pheno.ID =", Pheno.ID))
  # data = data.simu.surv.new(N = N, nSNP = nSNP, event.rate = event.rate,D:\R\1000Genome\Rcode\admixture-pop-structure6\data\phenotype GMat = Geno.mtx)
  data = data.simu.surv.admix(N=N, N1=N1, N2=N2, event.rate1 = ER1, event.rate2 = ER2, 
                              gammaVec = gammaVec, g.mtx = Geno.mtx.adj, g1.mtx = Geno.mtx.EUR.adj, g2.mtx = Geno.mtx.EAS.adj, a = a)
  
  Phen.mtx = data.frame(Pheno.ID = paste0("Pheno.ID-", Pheno.ID),
                        ER = paste0("ER in EUR = ", ER1, sep=", ", "ER in EAS = ", ER2),
                        ID = paste0("IID-",1:N),
                        event = data$event,
                        time = data$surv.time,
                        X1.new = data$X1.new,
                        X2 = data$X2,
                        X3 = data$X3)
  
  save(Phen.mtx, file = paste0("/gdata01/user/yuzhuoma/SPA-G/power/data/phenotype-homo-MAFsumover0.002-8log10MAF-neggamma/event.rate-",
                               ER, "/Phenotype-", Pheno.ID,".Rdata"))
  
} 
# }



