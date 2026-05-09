.libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))

library(doSNOW) 
#library(doParallel) 
#library(foreach) 
library(doRNG) 
library(brm) 
source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRR.R") 
source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRD.R") 
source("/home/yanjin41/brmplus_simulation/compare/1_CallMLE.R") 
source("/home/yanjin41/brmplus_simulation/compare/1.1_MLE_Point.R") 
source("/home/yanjin41/brmplus_simulation/compare/MLE_Point_Firth_for_RR.R") 
source("/home/yanjin41/brmplus_simulation/compare/MLE_Point_Firth_for_RD.R") 
source("/home/yanjin41/brmplus_simulation/compare/1.2_MLE_Var.R") 
source("/home/yanjin41/brmplus_simulation/compare/bayes_p.R") 
source("/home/yanjin41/brmplus_simulation/compare/MyFunc.R") 
source("/home/yanjin41/brmplus_simulation/compare/CI_exact_fast.R") 
source("/home/yanjin41/brmplus_simulation/compare/CI_LRT.R") 
source("/home/yanjin41/brmplus_simulation/compare/data_generation_simulation.R") 
#source("/home/yanjin41/brmplus_simulation/compare/data_ITE.R") 
source("/home/yanjin41/brmplus_simulation/R/RcppExports.R") 

library(epitools) 
library(geepack) 
library(sandwich) 
library(lmtest) 
library(brglm2) 
library(logistf) 
library(binom) 
library(epiR) 
library(PropCIs) 
library(MASS) 

### Modifiable parameters 
param = 'RR' # or 'RR' 
event = 'common' # 'rare' or 'common' 
hypothesis = 'null' # 'null' or 'alternative' 
argv <- commandArgs(TRUE) 
if (length(argv) == 0) { 
  print("No arguments supplied.") 
  n <- 50 
  R <- 200 
  } else { 
    for (i in 1:length(argv)) 
    eval(parse(text = argv[[i]])) 
    } 

#### 
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 
ncores <- min(ncores, 120) 
cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl) 

registerDoRNG(1234) 

result.mle <- foreach(r = (R-99):R, .packages = c("brm","epitools","geepack","sandwich","lmtest","brglm2", "MASS","logistf","binom","epiR","PropCIs"), .options.RNG=1234) %dopar% { 
        set.seed(r) 
        r1 <- run(param,n,event,hypothesis) 
        list(estimate = r1[1,],
             se = r1[2,], 
             low = r1[3,], 
             up = r1[4,], 
             p = r1[5,]) 
        } 

stopCluster(cl) 

result.all <- do.call(rbind, lapply(result.mle, as.data.frame)) 

write.csv(result.all, file = paste0(result_dir, "/np_simulation_results_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", R,".csv"))
