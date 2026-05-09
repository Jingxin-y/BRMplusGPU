.libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))



src_all <- function(){
  source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRR.R")
  source("/home/yanjin41/brmplus_simulation/compare/getProbScalarRD.R")
  source("/home/yanjin41/brmplus_simulation/compare/1_CallMLE.R")
  source("/home/yanjin41/brmplus_simulation/compare/1.1_MLE_Point.R")
  source("/home/yanjin41/brmplus_simulation/compare/1.2_MLE_Var.R")
  source("/home/yanjin41/brmplus_simulation/compare/bayes_p.R")
  source("/home/yanjin41/brmplus_simulation/compare/MyFunc.R")
  source("/home/yanjin41/brmplus_simulation/R/RcppExports.R")
  NULL
}
src_all()

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(doRNG)
  library(doSNOW)
  
library(brm)
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
})


param_vec <- c("RR")
event_vec <- c("rare")#"common",
hyp_vec   <- c("alternative")#，"null", 
ncores <- 120




argv <- commandArgs(TRUE) 
if (length(argv) == 0) { 
  print("No arguments supplied.") 
  n <- 50 
  R <- 200 
  } else { 
    for (i in 1:length(argv)) 
    eval(parse(text = argv[[i]])) 
    } 
    

firth_logbin_try <- function(dat,
                             start = rep(-0.01, 3),
                             eps = 1e-6,
                             maxit = 5000,
                             quiet = TRUE) {
  # return: c(est, se, lcl, ucl, p); NA if failure
  out_na <- function() setNames(rep(NA_real_, 5), c("est","se","lcl","ucl","p"))
  
  # basic sanity: required columns and finite values
  need <- c("y", "x", "v.1", "v.2")
  if (!all(need %in% names(dat))) return(out_na())
  dd <- dat[, need]
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (nrow(dd) == 0) return(out_na())
  for (nm in need) {
    if (any(!is.finite(dd[[nm]]))) return(out_na())
  }
  
  fit <- tryCatch(
    suppressWarnings(
      glm(y ~ x + v.1 + v.2 - 1,
          family  = binomial(link = "log"),
          data    = dd,
          start   = start,
          method  = "brglmFit",
          type    = "MPL_Jeffreys",
          control = brglmControl(epsilon = eps, maxit = maxit))
    ),
    error = function(e) {
      if (!quiet) message("firth_logbin_try error: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(fit)) return(out_na())
  
  # extract + guard against weird vcov / missing coefficient
  est <- tryCatch(unname(coef(fit)[["x"]]), error = function(e) NA_real_)
  V   <- tryCatch(vcov(fit), error = function(e) NULL)
  se  <- if (!is.null(V) && "x" %in% rownames(V) && "x" %in% colnames(V)) {
    sqrt(unname(V["x", "x"]))
  } else NA_real_
  
  if (!is.finite(est) || !is.finite(se) || se <= 0) return(out_na())
  
  z  <- est / se
  p  <- 2 * (1 - pnorm(abs(z)))
  ci <- est + c(-1.96, 1.96) * se
  
  c(est, se, ci[1], ci[2], p)
}

firth_logpois <- function(dat){
  fit <- glm(y ~ x + v.1 + v.2 - 1,
             family = poisson(link="log"),
             data   = dat,
             method = "brglmFit",
             type   = "MPL_Jeffreys",
             control = brglmControl(epsilon = 1e-6, maxit = 5000))
  
  est <- coef(fit)["x"]
  se  <- sqrt(vcov(fit)["x","x"])
  p   <- 2 * (1 - pnorm(abs(est/se)))
  ci  <- est + c(-1.96, 1.96) * se
  c(est, se, ci[1], ci[2], p)
}

firth_robust_logpois <- function(dat){
  fit <- glm(y ~ x + v.1 + v.2 - 1,
             family = poisson(link="log"),
             data   = dat,
             method = "brglmFit",
             type   = "MPL_Jeffreys",
             x      = TRUE,
             control = brglmControl(epsilon = 1e-6, maxit = 5000))
  
  est <- coef(fit)["x"]
  V   <- sandwich::vcovHC(fit, type="HC0")   # robust
  se  <- sqrt(V["x","x"])
  p   <- 2 * (1 - pnorm(abs(est/se)))
  ci  <- est + c(-1.96, 1.96) * se
  c(est, se, ci[1], ci[2], p)
}


firth_identity_rd_try <- function(dat,
                                  start = rep(-0.05, 3),
                                  eps_mu = 1e-8,
                                  maxit  = 5000,
                                  quiet  = TRUE) {
  # returns:
  # c(est, se, se_robust, lcl, ucl, lcl_robust, ucl_robust, p, p_robust)
  out_na <- function() {
    setNames(rep(NA_real_, 9),
             c("est","se","se_robust","lcl","ucl","lcl_robust","ucl_robust","p","p_robust"))
  }

  need <- c("y","x","v.1","v.2")
  if (!all(need %in% names(dat))) return(out_na())
  dd <- dat[, need, drop = FALSE]
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (nrow(dd) == 0) return(out_na())
  for (nm in need) if (any(!is.finite(dd[[nm]]))) return(out_na())

  X <- tryCatch(
    stats::model.matrix(stats::as.formula("~ x + v.1 + v.2 - 1"), data = dd),
    error = function(e) NULL
  )
  if (is.null(X) || any(!is.finite(X))) return(out_na())

  y <- dd$y
  if (any(!is.finite(y))) return(out_na())

  fit <- tryCatch(
    firth_binom_identity_fit(y = y, X = X, beta0 = start, eps_mu = eps_mu, maxit = maxit),
    error = function(e) {
      if (!quiet) message("firth_identity_rd_try error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit) || !isTRUE(fit$converged)) return(out_na())

  beta_hat <- fit$beta
  mu <- as.vector(X %*% beta_hat)

  # if mu out of bounds even at solution (can happen numerically), bail
  if (any(!is.finite(mu)) || any(mu <= eps_mu) || any(mu >= 1 - eps_mu)) return(out_na())

  # ---- SE (model-based): inverse observed information of *unpenalized* likelihood ----
  w <- 1 / (mu * (1 - mu))
  I <- crossprod(X, X * w)
  R <- tryCatch(chol(I), error = function(e) NULL)
  if (is.null(R)) return(out_na())
  V <- chol2inv(R)  # (X'WX)^{-1}

  # ---- Robust SE (sandwich) using individual score contributions ----
  # score_i = x_i * (y_i - mu_i)/(mu_i(1-mu_i))
  a <- (y - mu) / (mu * (1 - mu))          # n-vector
  S <- X * a                               # n x p matrix of score contributions
  meat <- crossprod(S)                     # sum s_i s_i'
  V_rob <- V %*% meat %*% V                # sandwich

  # pull "x" coefficient (fallback to col 1 if no name)
  j <- if ("x" %in% colnames(X)) match("x", colnames(X)) else 1L
  est <- beta_hat[j]

  se <- sqrt(V[j, j])
  se_rob <- sqrt(V_rob[j, j])

  if (!is.finite(est) || !is.finite(se) || se <= 0 ||
      !is.finite(se_rob) || se_rob <= 0) return(out_na())

  z  <- est / se
  zr <- est / se_rob

  p  <- 2 * (1 - stats::pnorm(abs(z)))
  pr <- 2 * (1 - stats::pnorm(abs(zr)))

  ci  <- est + c(-1.96, 1.96) * se
  cir <- est + c(-1.96, 1.96) * se_rob

  setNames(c(est, se, se_rob, ci[1], ci[2], cir[1], cir[2], p, pr),
           c("est","se","se_robust","lcl","ucl","lcl_robust","ucl_robust","p","p_robust"))
}
## =========================
## 1) Truth-setting helper
## =========================
simulate.rr <- function(n, event, hypothesis){
  
    if (event == "common"){
      if (hypothesis == "null"){
        alpha.true <- 0
        beta.true  <- c(1.5, 0.6)
        gamma.true <- c(0, 0)
      }else{
        alpha.true <- 0.3
        beta.true  <- c(1.65, 0.5)
        gamma.true <- c(0, 0)
      }
    }else{
      if (hypothesis == "null"){
        alpha.true <- 0
        beta.true  <- c(-4.7, 0.5)
        gamma.true <- c(0, 0)
      }else{
        alpha.true <- 0.7
        beta.true  <- c(-5.5, 0.5)
        gamma.true <- c(0, 0)
      }
    }
    
    data.simulation <- data.generation('RR', n, alpha.true, beta.true, gamma.true)
    
    va = as.matrix(data.simulation$data$v.1,ncol = 1)
    vb = cbind(data.simulation$data$v.1,data.simulation$data$v.2)
    y = data.simulation$data$y
    x = data.simulation$data$x
    Na0 = data.simulation$count[1]
    Na1 = data.simulation$count[2]
    N0_1 = data.simulation$count[3]
    N1_1 = data.simulation$count[4]
    P0 = N0_1/Na0
    P1 = N1_1/Na1
    
    pa = length(alpha.true)
    pb = length(beta.true)
    alpha.start = rep(0,pa)
    beta.start = rep(0,pb)
    
    weight = rep(1, length(y))
    max.step = min(pa * 20, 1000)
    thres = 1e-6
    thres.dicho = 1e-3
    
    ##brm
    est.brm <- MLEst('RR', y, x, va, vb, weight, max.step, thres, alpha.start = rep(0, pa),
                     beta.start = rep(0, pb), pa, pb)
    
    ##CMH
    sam.CMH <- matrix(c(Na0-N0_1,Na1-N1_1,N0_1,N1_1),2,2)
    est.CMH <- riskratio(sam.CMH, method="small", correction=TRUE)
    
    ##
    v.1 = vb[,1]
    v.2 = vb[,2]
    ##log-binomial
    # est.lb1 <- glm(y~x+v.1+v.2-1, family = binomial(link = "log"), data = data.simulation$data, start = rep(-0.01,3))
    est.lb <- firth_logbin_try(data.simulation$data)
    
    
    ##log-poisson
    # est.lp1 <- glm(y~x+v.1+v.2-1, family = poisson(link = "log"), data = data.simulation$data)
    est.lp <- firth_logpois(data.simulation$data)
    
    ##robust log-poisson
    # est.rlp1 <- quasi.poisson(data.simulation$data)
    est.rlp <- firth_robust_logpois(data.simulation$data)
    
  
  
  ##brm + firth
  
  # est.brm.firth <- MLEst('RR', y, x, va, vb, weight, max.step, thres, alpha.start = rep(0, pa),
  #                        beta.start = rep(0, pb), pa, pb, method="firth")
  # 
  ## brm_ad
  est.brm.ad = est.brm
  if(P0==0|P0==1|P1==0|P1==1) {
    est.bayes = bayes_est_RR(Na0,Na1,N0_1,N1_1)
    est.brm.ad$point.est[1] = est.bayes$point.est
    est.brm.ad$se.est[1] = est.bayes$se.est
    est.brm.ad$conf.lower[1] = est.bayes$conf.lower
    est.brm.ad$conf.upper[1] = est.bayes$conf.upper
    est.brm.ad$p.value[1] = est.bayes$p.value
  }
  
  ###result
  point.est <- as.vector(c(est.brm$point.est[1],
                           est.brm.ad$point.est[1],
                           log(est.CMH$measure[2,1]),
                           est.lb[1],
                           est.lp[1],
                           est.rlp[1],
                           1,
                           1,
                           1,
                           1,
                           1,
                           1,
                           1,
                           1))
  se.est <- as.vector(c(est.brm$se.est[1],
                        est.brm.ad$se.est[1],
                        (log(est.CMH$measure[2,1])-log(est.CMH$measure[2,2]))/qnorm(0.975),
                        est.lb[2],
                        est.lp[2],
                        est.rlp[2],
                        1,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1))
  con.lower <- as.vector(c(est.brm$conf.lower[1],
                           est.brm.ad$conf.lower[1],
                           log(est.CMH$measure[2,2]),
                           est.lb[3],
                           est.lp[3],
                           est.rlp[3],
                           1,
                           1,
                           1,
                           1,
                           1,
                           1,
                           1,
                           1))
  con.upper <- as.vector(c(est.brm$conf.upper[1],
                           est.brm.ad$conf.upper[1],
                           log(est.CMH$measure[2,3]),
                           est.lb[4],
                           est.lp[4],
                           est.rlp[4],
                           1,
                           1,
                           1,
                           1,
                           1,
                           1,
                           1,
                           1))
  
  p.value <- as.vector(c(est.brm$p.value[1],
                         est.brm.ad$p.value[1],
                         est.CMH$p.value[2,1],
                         est.lb[5],
                         est.lp[5],
                         est.rlp[5],
                         1,
                         1,
                         1,
                         1,
                         1,
                         1,
                         1,
                         1))
  
  result.comp <- rbind(point.est,se.est,con.lower,con.upper,p.value)
  colnames(result.comp) <- c("brm","brm_ad","CMH","log-binomial","log-poisson","robust log-possion","brm_firth",
                             "brm_exact","brm_exact_ad","g-computation","GC_BR","GC_FC","GC_FC_BR1","GC_FC_BR2")
  return(result.comp)
}

simulate.rd <- function(n, event, hypothesis){
  
    if (event == "common"){
      if (hypothesis == "null"){
        alpha.true = 0
        beta.true   = c(0.9,0.5)
        gamma.true  = c(0,0)
      }else{
        alpha.true = 0.1
        beta.true   = c(0.9,0.2)
        gamma.true  = c(0,0)
      }
    }else{
      if (hypothesis == "null"){
        alpha.true = 0
        beta.true   = c(-4.5,0.5)
        gamma.true  = c(0,0)
      }else{
        alpha.true = 0.05
        beta.true   = c(-5.5,0.2)
        gamma.true  = c(0,0)# rare
      }
    }
    
    data.simulation <- data.generation('RD', n, alpha.true, beta.true, gamma.true)
    va = as.matrix(data.simulation$data$v.1,ncol = 1)
    vb = cbind(data.simulation$data$v.1,data.simulation$data$v.2)
    y = data.simulation$data$y
    x = data.simulation$data$x
    Na0 = data.simulation$count[1]
    Na1 = data.simulation$count[2]
    N0_1 = data.simulation$count[3]
    N1_1 = data.simulation$count[4]
    
    P0 = N0_1/Na0
    P1 = N1_1/Na1
    
    pa = length(alpha.true)
    pb = length(beta.true)
    alpha.start = rep(0,pa)
    beta.start = rep(0,pb)
    
    weight = rep(1, length(y))
    max.step = min(pa * 20, 1000)
    thres = 1e-6
    
    ##brm
    est.brm.or <- MLEst('RD', y, x, va, vb, weight, max.step, thres, alpha.start = rep(0, pa),
                        beta.start = rep(0, pb), pa, pb)
    
    est.brm.ad = est.brm.or
    if(P0==0|P0==1|P1==0|P1==1) {
      est.bayes = bayes_est_RD(Na0,Na1,N0_1,N1_1)
      est.brm.ad$point.est[1] = est.bayes$point.est
      est.brm.ad$se.est[1] = est.bayes$se.est
      est.brm.ad$conf.lower[1] = est.bayes$conf.lower
      est.brm.ad$conf.upper[1] = est.bayes$conf.upper
      est.brm.ad$p.value[1] = est.bayes$p.value
    }
    
    ## bayesian prior
    est.bayes = bayes_est_RD(Na0,Na1,N0_1,N1_1)
    
    v.1 = vb[,1]
    v.2 = vb[,2]
    ## GLM with identity link (calc_risk with identity link?)
    e.glm <- firth_identity_rd_try(data.simulation$data)
  
  ###result
  point.est <- c(est.brm.or$point.est[1],
                 est.brm.ad$point.est[1],
                 est.bayes$point.est,
                 e.glm["est"],
                 e.glm["est"],
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1)
  se.est <- c(est.brm.or$se.est[1],
              est.brm.ad$se.est[1],
              est.bayes$se.est,
              e.glm["se"],
              e.glm["se_robust"],
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              1)
  CI.low.or <- c(est.brm.or$conf.lower[1],
                 est.brm.ad$conf.lower[1],
                 est.bayes$conf.lower,
                 e.glm["lcl"],
                 e.glm["lcl_robust"],
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1)
  CI.up.or <- c(est.brm.or$conf.upper[1],
                est.brm.ad$conf.upper[1],
                est.bayes$conf.upper,
                e.glm["ucl"],
                e.glm["ucl_robust"],
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1)
  p.value <- c(est.brm.or$p.value[1],
               est.brm.ad$p.value[1],
               est.bayes$p.value,
               e.glm["p"],
               e.glm["p_robust"],
               1,
               1,
               1,
               1,
               1,
               1,
               1,
               1,
               1)
  
  result.comp <- rbind(point.est,se.est,CI.low.or,CI.up.or,p.value)
  colnames(result.comp) <- c("brm","brm_ad","bayes","glm","lpm","MN", "firth",
                             "brm_exact","brm_exact_ad","g-computation","GC_BR","GC_FC","GC_FC_BR1","GC_FC_BR2")
  return(result.comp)
}


## =========================
## 2) One replicate (returns named pvals)
## =========================
one_rep <- function(r, param, n, event, hypothesis,
                    max.step = NULL, thres = 1e-6){

  set.seed(r)

  simulate.fun <- if (param == "RR") simulate.rr else simulate.rd

  full <- tryCatch(
    simulate.fun(n, event, hypothesis),
    error = function(e) {
      # 出错就返回一个全 NA 的 5xK 占位矩阵（K 只能用已知列名来定）
      # 这里最稳的做法：直接返回 NULL，让上层处理
      NULL
    }
  )

  ok <- !is.null(full) && any(is.finite(full))

  list(ok = ok, full = full)
}

## =========================
## 3) Run one scenario -> RETURN p_mat + (optional) SAVE p_mat
## =========================
run_scenario <- function(param, n, event, hypothesis, R,
                         ncores = 8, thres = 1e-6, max.step = NULL,
                         result_dir = "/scratch/yanjin41/RRRDOR/brmplus_simulation", save_mat = TRUE){
  
  ## cluster
  cl <- makeCluster(ncores)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)
  registerDoParallel(cl)
  
  ## 让每个 worker 都有同样的库路径 + 包 + source
  clusterEvalQ(cl, {
    .libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))
    suppressPackageStartupMessages({
     library(brm);library(epitools);library(geepack);library(sandwich)
     library(lmtest);library(brglm2);library(logistf);library(binom)
     library(epiR);library(PropCIs);library(MASS)
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
    r = (R-999):R,   # 正式跑用 1:R
    .export = c("one_rep", "firth_logbin_try", "firth_logpois",
              "firth_robust_logpois", "firth_identity_rd_try",
              "simulate.rr", "simulate.rd"),
    .noexport = c()
  ) %dopar% {
    one_rep(r, param, n, event, hypothesis, max.step = max.step, thres = thres)
  }

  ok_vec <- vapply(res, `[[`, logical(1), "ok")

  # --- 把每次 replicate 的 full (5 x K) 堆叠成 full_arr: 5 x K x Rrun ---
  full_list <- lapply(res, `[[`, "full")

  # 假设每次 full 的维度一致（5 x K）
  d1 <- nrow(full_list[[1]])
  d2 <- ncol(full_list[[1]])
  Rrun <- length(full_list)

  full_arr <- array(NA_real_, dim = c(d1, d2, Rrun))
  for(i in seq_len(Rrun)){
    full_arr[,,i] <- full_list[[i]]
  }
  dimnames(full_arr) <- list(
    rownames(full_list[[1]]),
    colnames(full_list[[1]]),
    paste0("rep_", seq_len(Rrun))
  )

  if (save_mat) {
    if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    tag <- paste0("param=",param,
                "_n=",n,
                "_event=",event,
                "_hyp=",hypothesis,
                "_R=",R)

  # 2) 保存 full_arr（保存所有输出）
    saveRDS(full_arr, file = file.path(result_dir, paste0("full_arr_", tag, ".rds")))

  }
  list(full_arr = full_arr, ok = ok_vec)
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
    save_mat = TRUE
  )
}
print(Sys.time() - t0)
