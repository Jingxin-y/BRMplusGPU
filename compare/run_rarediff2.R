.libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))



src_all <- function(){
  source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRR.R")
  source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRD.R")
  source("/home/yanjin41/brmplus_simulation/compare/1_CallMLE.R")
  source("/home/yanjin41/brmplus_simulation/compare/1.1_MLE_Point.R")
  source("/home/yanjin41/brmplus_simulation/compare/1.2_MLE_Var.R")
  source("/home/yanjin41/brmplus_simulation/compare/bayes_p.R")
  source("/home/yanjin41/brmplus_simulation/compare/MyFunc.R")
  source("/home/yanjin41/brmplus_simulation/compare/CI_exact_diff.R")
  source("/home/yanjin41/brmplus_simulation/compare/data_generation_simulation.R")
  source("/home/yanjin41/brmplus_simulation/R/RcppExports.R")
  NULL
}
src_all()

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(doRNG)
  library(doSNOW)
  library(lmtest)
  
  library(PropCIs)
  library(epitools)
  library(brm)
  library(MASS)
  library(sandwich)
})


param_vec <- c("RR", "RD")
event_vec <- "rare11"
hyp_vec   <- "alternative"
ncores <- 100

argv <- commandArgs(TRUE) 
if (length(argv) == 0) { 
  print("No arguments supplied.") 
  n <- 50 
  R <- 200 
  } else { 
    for (i in 1:length(argv)) 
    eval(parse(text = argv[[i]])) 
    } 



## =========================
## 1) Truth-setting helper
## =========================
get_truth <- function(param, event, hypothesis){
  if(param == "RR"){
        alpha.true <- 2.4
        beta.true  <- c(-6.54, -0.5)
        gamma.true <- c(0.2, -0.5)
  } else { # RD
        alpha.true <- 0.1
        beta.true  <- c(-6.54, -0.5)
        gamma.true <- c(0.2, -0.5)
  }
  list(alpha.true = alpha.true, beta.true = beta.true, gamma.true = gamma.true)
}

## =========================
## 2) One replicate (returns named pvals)
## =========================
one_rep <- function(r, param, n, event, hypothesis,
                    max.step = NULL, thres = 1e-6){
  
  set.seed(r)
  
  tru <- get_truth(param, event, hypothesis)
  alpha.true <- tru$alpha.true
  beta.true  <- tru$beta.true
  gamma.true <- tru$gamma.true
  
  dat <- data.generation(param, n, alpha.true, beta.true, gamma.true)
  
  y  <- dat$data$y
  x  <- dat$data$x
  va <- as.matrix(dat$data$v.1, ncol = 1)
  vb <- cbind(dat$data$v.1, dat$data$v.2)
  
  Na0  <- dat$count[1]
  Na1  <- dat$count[2]
  N0_1 <- dat$count[3]
  N1_1 <- dat$count[4]
  
  pa <- length(alpha.true)
  pb <- length(beta.true)
  
  alpha.start <- rep(0, pa)
  beta.start  <- rep(0, pb)
  
  weight <- rep(1, length(y))
  if (is.null(max.step)) max.step <- min(pa * 20, 1000)
  
  out <- tryCatch({
  
    P0 <- N0_1 / Na0
    P1 <- N1_1 / Na1
    
    ## brm MLE + BC exact
    fit_brm_RR <- MLEst("RR", y, x, va, vb, weight, max.step, thres,
                        alpha.start = alpha.start, beta.start = beta.start, pa, pb)
                        
    RR.brm.ad = fit_brm_RR
  if(P0==0|P0==1|P1==0|P1==1) {
    est.bayes = bayes_est_RR(Na0,Na1,N0_1,N1_1)
    RR.brm.ad$point.est[1] = est.bayes$point.est
    RR.brm.ad$se.est[1] = est.bayes$se.est
    RR.brm.ad$conf.lower[1] = est.bayes$conf.lower
    RR.brm.ad$conf.upper[1] = est.bayes$conf.upper
    RR.brm.ad$p.value[1] = est.bayes$p.value
  }
    exact_RR <- exact("RR", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      RR.brm.ad$point.est, RR.brm.ad$se.est, pa, pb)
    
    fit_brm_RD <- MLEst("RD", y, x, va, vb, weight, max.step, thres,
                        alpha.start = alpha.start, beta.start = beta.start, pa, pb)
                        
     RD.brm.ad = fit_brm_RD
  if(P0==0|P0==1|P1==0|P1==1) {
    est.bayes = bayes_est_RD(Na0,Na1,N0_1,N1_1)
    RD.brm.ad$point.est[1] = est.bayes$point.est
    RD.brm.ad$se.est[1] = est.bayes$se.est
    RD.brm.ad$conf.lower[1] = est.bayes$conf.lower
    RD.brm.ad$conf.upper[1] = est.bayes$conf.upper
    RD.brm.ad$p.value[1] = est.bayes$p.value
  }
    exact_RD <- exact("RD", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      RD.brm.ad$point.est, RD.brm.ad$se.est, pa, pb)
    
    ## regressions
    v.1 <- vb[, 1]
    v.2 <- vb[, 2]
    
    fit_lb <- glm(y ~ x + v.1 + v.2 - 1,
                  family = binomial(link = "log"),
                  data   = dat$data,
                  start  = rep(-0.01, 3))
    
    fit_lp <- glm(y ~ x + v.1 + v.2 - 1,
                  family = poisson(link = "log"),
                  data   = dat$data)
    
    fit_rlp <- quasi.poisson(dat$data)
    
    fit_glm_id <- glm(y ~ x + v.1 + v.2 - 1,
                      family = binomial(link = "identity"),
                      data   = dat$data,
                      start  = rep(0.01, 3))
    
     fit_lpm <- lm(y ~ x + v.1 + v.2 - 1, data = dat$data)
     e.lpm <- coeftest(fit_lpm, vcov = vcovHC(fit_lpm, type = "HC3"))

    
    p.lpm <- e.lpm[1, 4]
    
    ## MN + CMH (table-based)
    
    
    p_mn <- tryCatch({
      mn_point <- P1 - P0
      mn_ci <- PropCIs::diffscoreci(N1_1, Na1, N0_1, Na0, conf.level = 0.95)
      mn_se <- (mn_ci$conf.int[2] - mn_ci$conf.int[1]) / (2 * qnorm(0.975))
      if (!is.finite(mn_se) || mn_se <= 0) return(NA_real_)
      z <- mn_point / mn_se
      2 * pmin(pnorm(z), 1 - pnorm(z))
    }, error = function(e) NA_real_)
    
    p_cmh <- tryCatch({
      sam_2x2 <- matrix(c(Na0 - N0_1, Na1 - N1_1,
                          N0_1,       N1_1),
                        nrow = 2, byrow = FALSE)
      est_cmh <- epitools::riskratio(sam_2x2, method = "small", correction = TRUE)
      as.numeric(est_cmh$p.value[2, 1])
    }, error = function(e) NA_real_)
    
    pvals <- c(
      brm_RR    = fit_brm_RR$p.value[1],
      brm_RD    = fit_brm_RD$p.value[1],
      brm_BC_RR = exact_RR$p[1],
      brm_BC_RD = exact_RD$p[1],
      lb        = summary(fit_lb)$coefficients[1, 4],
      lp        = summary(fit_lp)$coefficients[1, 4],
      rlp       = as.numeric(fit_rlp[5]),
      glm_id    = summary(fit_glm_id)$coefficients[1, 4],
      lpm       = p.lpm,
      MN        = p_mn,
      CMH       = p_cmh
    )
    
    list(ok = TRUE, p = pvals, err = NA_character_)
    
  }, error = function(e){
    list(
      ok = FALSE,
      p  = c(brm_RR=NA_real_, brm_RD=NA_real_, brm_BC_RR=NA_real_, brm_BC_RD=NA_real_,
             lb=NA_real_, lp=NA_real_, rlp=NA_real_, glm_id=NA_real_, lpm=NA_real_,
             MN=NA_real_, CMH=NA_real_),
      err = conditionMessage(e)
    )
  })
  
  out
}


## =========================
## 3) Run one scenario -> RETURN p_mat + (optional) SAVE p_mat
## =========================
run_scenario <- function(param, n, event, hypothesis, R,
                         ncores = 8, thres = 1e-6, max.step = NULL,
                         result_dir = "/scratch/yanjin41/RRRDOR/brmplus_simulation", save_pmat = TRUE){
  
  ## cluster
  cl <- makeCluster(ncores)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)
  registerDoParallel(cl)
  
  ## 让每个 worker 都有同样的库路径 + 包 + source
  clusterEvalQ(cl, {
    .libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))
    suppressPackageStartupMessages({
      library(PropCIs); library(epitools); library(brm)
      library(MASS); library(sandwich); library(lmtest)
      library(doRNG); library(foreach)
    })
    NULL
  })
  

  
  clusterEvalQ(cl, {
    source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRR.R")
    source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRD.R")
    source("/home/yanjin41/brmplus_simulation/compare/1_CallMLE.R")
    source("/home/yanjin41/brmplus_simulation/compare/1.1_MLE_Point.R")
    source("/home/yanjin41/brmplus_simulation/compare/1.2_MLE_Var.R")
    source("/home/yanjin41/brmplus_simulation/compare/bayes_p.R")
    source("/home/yanjin41/brmplus_simulation/compare/MyFunc.R")
    source("/home/yanjin41/brmplus_simulation/compare/CI_exact_diff.R")
    source("/home/yanjin41/brmplus_simulation/compare/data_generation_simulation.R")
    source("/home/yanjin41/brmplus_simulation/R/RcppExports.R")
    NULL
  })
  
  ## 随机数：稳定复现
  doRNG::registerDoRNG(1234)
  
  res <- foreach(
    r = (R-999):R,
    .export = c("one_rep", "get_truth"),
    .noexport = c()   # 明确不需要其它对象
  ) %dopar% {
    one_rep(r, param, n, event, hypothesis, max.step = max.step, thres = thres)
  }
  
  ok_vec <- vapply(res, `[[`, logical(1), "ok")
  
  p_mat <- do.call(rbind, lapply(res, `[[`, "p"))
  colnames(p_mat) <- c("brm_RR","brm_RD","brm_BC_RR","brm_BC_RD",
                       "lb","lp","rlp","glm_id","lpm","MN","CMH")
  
  meta <- data.frame(
    n = n, event = event, hypothesis = hypothesis, param = param,
    R = R,
    success_rate = mean(ok_vec),
    stringsAsFactors = FALSE
  )
  
  ## 保存 p_mat（每个 scenario 一个文件）
  if (save_pmat) {
    if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    tag <- paste0("pmat_param=",param,
                  "_n=",n,
                  "_event=",event,
                  "_hyp=",hypothesis,
                  "_R=",R)
    write.csv(p_mat, file = paste0(result_dir, "/results_",tag,".csv"))

  }
  
  ## 同时返回给你（交互跑也可直接拿到）
  list(meta = meta, p_mat = p_mat, ok = ok_vec)
}




scenarios <- expand.grid(
  n = n,
  event = event_vec,
  hypothesis = hyp_vec,
  param = param_vec,
  stringsAsFactors = FALSE
)


all_out <- vector("list", nrow(scenarios))

t0 <- Sys.time()
for(s in seq_len(nrow(scenarios))){
  cat("Running scenario", s, "of", nrow(scenarios), "...\n")
  all_out[[s]] <- run_scenario(
    param = scenarios$param[s],
    n = scenarios$n[s],
    event = scenarios$event[s],
    hypothesis = scenarios$hypothesis[s],
    R = R,
    ncores = ncores,
    save_pmat = TRUE
  )
}
print(Sys.time() - t0)
