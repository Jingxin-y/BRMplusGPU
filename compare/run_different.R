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
  
  library(PropCIs)
  library(epitools)
  library(brm)
  library(MASS)
  library(sandwich)
  library(geepack)
  library(lmtest)
  library(brglm2)
  library(logistf)
  library(binom)
  library(epiR)
  library(PropCIs)
})


param_vec <- c("RR", "RD")
event_vec <- c("rare")
hyp_vec   <- c("null", "alternative")
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
    if (event == "common") {
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true  <- c(1.5, 0.6)
        gamma.true <- c(0.2, -0.5)
      } else {
        alpha.true <- 0.3
        beta.true  <- c(1.65, 0.5)
        gamma.true <- c(0.2, -0.5)
      }
    } else { # rare
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true  <- c(-4.7, 0.5)
        gamma.true <- c(0.2, -0.5)
      } else {
        alpha.true <- 0.7
        beta.true  <- c(-5.5, 0.5)
        gamma.true <- c(0.2, -0.5)
      }
    }
  } else { # RD
    if (event == "common") {
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true  <- c(0.9, 0.5)
        gamma.true <- c(0.2, -0.5)
      } else {
        alpha.true <- 0.1
        beta.true  <- c(0.9, 0.2)
        gamma.true <- c(0.2, -0.5)
      }
    } else { # rare
      if (hypothesis == "null") {
        alpha.true <- 0
        beta.true  <- c(-4.5, 0.5)
        gamma.true <- c(0.2, -0.5)
      } else {
        alpha.true <- 0.05
        beta.true  <- c(-5.5, 0.2)
        gamma.true <- c(0.2, -0.5)
      }
    }
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
    
    ## brm MLE + BC exact
    fit_brm_RR <- MLEst("RR", y, x, va, vb, weight, max.step, thres,
                        alpha.start = alpha.start, beta.start = beta.start, pa, pb)
    exact_RR <- exact("RR", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      fit_brm_RR$point.est, fit_brm_RR$se.est, pa, pb)
    
    fit_brm_RD <- MLEst("RD", y, x, va, vb, weight, max.step, thres,
                        alpha.start = alpha.start, beta.start = beta.start, pa, pb)
    exact_RD <- exact("RD", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      fit_brm_RD$point.est, fit_brm_RD$se.est, pa, pb)
    
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
    
    ## MN + CMH (table-based)
    P0 <- N0_1 / Na0
    P1 <- N1_1 / Na1
    
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
    
    
    ## g-computaion & g-computation_BR
Y1 <- y[which(x == 1)]
Y0 <- y[which(x == 0)]
V2.1 <- v.2[which(x == 1)]
V2.0 <- v.2[which(x == 0)]
X1 <- x[which(x == 1)]
X0 <- x[which(x == 0)]

data.treat <- data.frame(Y1, V2.1)
data.control <- data.frame(Y0, V2.0)

est.treat <- glm(Y1 ~ V2.1, family = binomial, data = data.treat)
est.control <- glm(Y0 ~ V2.0, family = binomial, data = data.control)

beta.hat.treat <- est.treat$coefficients
beta.hat.control <- est.control$coefficients

V.FC.treat <- cbind(1, V2.1)
V.FC.control <- cbind(1, V2.0)


beta.hat.star.treat <- beta.hat.treat + colMeans(hatvalues(est.treat) * phi(Y1, V.FC.treat, beta.hat.treat, sum(x == 1) / n))
beta.hat.star.control <- beta.hat.control + colMeans(hatvalues(est.control) * phi(Y0, V.FC.control, beta.hat.control, sum(x == 0) / n))

# beta_hat
p.hat.treat <- mean(c(Y1, predict(est.treat, newdata = data.control, type = "response")))
p.hat.control <- mean(c(Y0, predict(est.control, newdata = data.treat, type = "response")))

# beta_hat_star
p.hat.star.treat <- mean(c(Y1, m(V.FC.control %*% beta.hat.star.treat)))
p.hat.star.control <- mean(c(Y0, m(V.FC.treat %*% beta.hat.star.control)))

# beta_tilde
fit.treat <- logistf(Y1 ~ V2.1, data = data.treat)
fit.control <- logistf(Y0 ~ V2.0, data = data.control)

beta.tilde.treat <- fit.treat$coefficients
beta.tilde.control <- fit.control$coefficients

# beta_tilde_star
beta.tilde.star.treat <- beta.tilde.treat + colMeans(as.vector(hii(V.FC.treat, beta.tilde.treat)) * (phi(Y1, V.FC.treat, beta.tilde.treat, sum(x == 1) / n)
                                                                                                     - (V.FC.treat * as.vector(1 - 2 * m(V.FC.treat %*% beta.tilde.treat))) %*% t(ginv(fish(V.FC.treat, beta.tilde.treat))) / 2))
beta.tilde.star.control <- beta.tilde.control + colMeans(as.vector(hii(V.FC.control, beta.tilde.control)) * (phi(Y0, V.FC.control, beta.tilde.control, sum(x == 0) / n)
                                                                                                             - (V.FC.control * as.vector(1 - 2 * m(V.FC.control %*% beta.tilde.control))) %*% t(ginv(fish(V.FC.control, beta.tilde.control))) / 2))
# beta_tilde_doustar
beta.tilde.doustar.treat <- beta.tilde.treat - colMeans(as.vector(hii(V.FC.treat, beta.tilde.treat)) * ((V.FC.treat * as.vector(1 - 2 * m(V.FC.treat %*% beta.tilde.treat))) %*% t(ginv(fish(V.FC.treat, beta.tilde.treat))) / 2))
beta.tilde.doustar.control <- beta.tilde.control - colMeans(as.vector(hii(V.FC.control, beta.tilde.control)) * ((V.FC.control * as.vector(1 - 2 * m(V.FC.control %*% beta.tilde.control))) %*% t(ginv(fish(V.FC.control, beta.tilde.control))) / 2))

# beta_tilde
p.tilde.treat <- mean(c(Y1, m(V.FC.control %*% beta.tilde.treat)))
p.tilde.control <- mean(c(Y0, m(V.FC.treat %*% beta.tilde.control)))

# beta_tilde_star
p.tilde.star.treat <- mean(c(Y1, m(V.FC.control %*% beta.tilde.star.treat)))
p.tilde.star.control <- mean(c(Y0, m(V.FC.treat %*% beta.tilde.star.control)))

# beta_tilde_starstar
p.tilde.doustar.treat <- mean(c(Y1, m(V.FC.control %*% beta.tilde.doustar.treat)))
p.tilde.doustar.control <- mean(c(Y0, m(V.FC.treat %*% beta.tilde.doustar.control)))

li.hat <- l.mu(Y1, V.FC.treat, beta.hat.treat, Y0, V.FC.control, beta.hat.control)
li.hat.star <- l.mu(Y1, V.FC.treat, beta.hat.star.treat, Y0, V.FC.control, beta.hat.star.control)
li.tilde <- l.mu(Y1, V.FC.treat, beta.tilde.treat, Y0, V.FC.control, beta.tilde.control)
li.tilde.star <- l.mu(Y1, V.FC.treat, beta.tilde.star.treat, Y0, V.FC.control, beta.tilde.star.control)
li.tilde.doustar <- l.mu(Y1, V.FC.treat, beta.tilde.doustar.treat, Y0, V.FC.control, beta.tilde.doustar.control)

alpha.hat.RR <- log(p.hat.treat / p.hat.control)
alpha.hat.RD <- p.hat.treat - p.hat.control
alpha.hat.star.RR <- log(p.hat.star.treat / p.hat.star.control)
alpha.hat.star.RD <- p.hat.star.treat - p.hat.star.control
alpha.tilde.RR <- log(p.tilde.treat / p.tilde.control)
alpha.tilde.RD <- p.tilde.treat - p.tilde.control
alpha.tilde.star.RR <- log(p.tilde.star.treat / p.tilde.star.control)
alpha.tilde.star.RD <- p.tilde.star.treat - p.tilde.star.control
alpha.tilde.doustar.RR <- log(p.tilde.doustar.treat / p.tilde.doustar.control)
alpha.tilde.doustar.RD <- p.tilde.doustar.treat - p.tilde.doustar.control



se.hat.RR <- sqrt(var.est.RR(li.hat, p.hat.control, p.hat.treat))
se.hat.star.RR <- sqrt(var.est.RR(li.hat.star, p.hat.star.control, p.hat.star.treat))
se.tilde.RR <- sqrt(var.est.RR(li.tilde, p.tilde.control, p.tilde.treat))
se.tilde.star.RR <- sqrt(var.est.RR(li.tilde.star, p.tilde.star.control, p.tilde.star.treat))
se.tilde.doustar.RR <- sqrt(var.est.RR(li.tilde.doustar, p.tilde.doustar.control, p.tilde.doustar.treat))

se.hat.RD <- sqrt(var.est.RD(li.hat, p.hat.control, p.hat.treat))
se.hat.star.RD <- sqrt(var.est.RD(li.hat.star, p.hat.star.control, p.hat.star.treat))
se.tilde.RD <- sqrt(var.est.RD(li.tilde, p.tilde.control, p.tilde.treat))
se.tilde.star.RD <- sqrt(var.est.RD(li.tilde.star, p.tilde.star.control, p.tilde.star.treat))
se.tilde.doustar.RD <- sqrt(var.est.RD(li.tilde.doustar, p.tilde.doustar.control, p.tilde.doustar.treat))


p.hat.RR <- 2 * min(pnorm(alpha.hat.RR / se.hat.RR), 1 - pnorm(alpha.hat.RR / se.hat.RR))
p.hat.star.RR <- 2 * min(pnorm(alpha.hat.star.RR / se.hat.star.RR), 1 - pnorm(alpha.hat.star.RR / se.hat.star.RR))
p.tilde.RR <- 2 * min(pnorm(alpha.tilde.RR / se.tilde.RR), 1 - pnorm(alpha.tilde.RR / se.tilde.RR))
p.tilde.star.RR <- 2 * min(pnorm(alpha.tilde.star.RR / se.tilde.star.RR), 1 - pnorm(alpha.tilde.star.RR / se.tilde.star.RR))
p.tilde.doustar.RR <- 2 * min(pnorm(alpha.tilde.doustar.RR / se.tilde.doustar.RR), 1 - pnorm(alpha.tilde.doustar.RR / se.tilde.doustar.RR))

p.hat.RD <- 2 * min(pnorm(alpha.hat.RD / se.hat.RD), 1 - pnorm(alpha.hat.RD / se.hat.RD))
p.hat.star.RD <- 2 * min(pnorm(alpha.hat.star.RD / se.hat.star.RD), 1 - pnorm(alpha.hat.star.RD / se.hat.star.RD))
p.tilde.RD <- 2 * min(pnorm(alpha.tilde.RD / se.tilde.RD), 1 - pnorm(alpha.tilde.RD / se.tilde.RD))
p.tilde.star.RD <- 2 * min(pnorm(alpha.tilde.star.RD / se.tilde.star.RD), 1 - pnorm(alpha.tilde.star.RD / se.tilde.star.RD))
p.tilde.doustar.RD <- 2 * min(pnorm(alpha.tilde.doustar.RD / se.tilde.doustar.RD), 1 - pnorm(alpha.tilde.doustar.RD / se.tilde.doustar.RD))

    
    pvals <- c(
      brm_RR    = fit_brm_RR$p.value[1],
      brm_RD    = fit_brm_RD$p.value[1],
      brm_BC_RR = exact_RR$p[1],
      brm_BC_RD = exact_RD$p[1],
      lb        = summary(fit_lb)$coefficients[1, 4],
      lp        = summary(fit_lp)$coefficients[1, 4],
      rlp       = as.numeric(fit_rlp[5]),
      glm_id    = summary(fit_glm_id)$coefficients[1, 4],
      lpm       = summary(fit_lpm)$coefficients[1, 4],
      MN        = p_mn,
      CMH       = p_cmh,
      p.hat.RR = p.hat.RR,
      p.hat.star.RR = p.hat.star.RR,
      p.tilde.RR = p.tilde.RR,
      p.tilde.star.RR = p.tilde.star.RR,
      p.tilde.doustar.RR = p.tilde.doustar.RR,
      p.hat.RD = p.hat.RD,
      p.hat.star.RD = p.hat.star.RD,
      p.tilde.RD = p.tilde.RD,
      p.tilde.star.RD = p.tilde.star.RD,
      p.tilde.doustar.RD = p.tilde.doustar.RD
    )
    
    list(ok = TRUE, p = pvals, err = NA_character_)
    
  }, error = function(e){
  list(
    ok = FALSE,
    p  = c(brm_RR = NA_real_,brm_RD = NA_real_,brm_BC_RR = NA_real_,brm_BC_RD = NA_real_,
      lb = NA_real_,lp = NA_real_,rlp = NA_real_,glm_id = NA_real_,lpm = NA_real_,
      MN = NA_real_,CMH = NA_real_,p.hat.RR = NA_real_,p.hat.star.RR = NA_real_,
      p.tilde.RR = NA_real_,p.tilde.star.RR = NA_real_, p.tilde.doustar.RR = NA_real_,
      p.hat.RD = NA_real_,p.hat.star.RD = NA_real_,p.tilde.RD = NA_real_,
      p.tilde.star.RD = NA_real_,p.tilde.doustar.RD = NA_real_
    ),
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
      library(MASS); library(sandwich);library(geepack)
      library(lmtest);library(brglm2);library(logistf)
      library(binom);library(epiR)
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
    r = (R-499):R,
    .export = c("one_rep", "get_truth"),
    .noexport = c()   # 明确不需要其它对象
  ) %dopar% {
    one_rep(r, param, n, event, hypothesis, max.step = max.step, thres = thres)
  }
  
  ok_vec <- vapply(res, `[[`, logical(1), "ok")
  
  p_mat <- do.call(rbind, lapply(res, `[[`, "p"))
  colnames(p_mat) <- c("brm_RR","brm_RD","brm_BC_RR","brm_BC_RD",
                       "lb","lp","rlp","glm_id","lpm","MN","CMH",
                       "GC_RR","GCBR_RR","GCFC_RR","GCFCBR1_RR","GCFCBR2_RR",
                       "GC_RD","GCBR_RD","GCFC_RD","GCFCBR1_RD","GCFCBR2_RD")
  
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
