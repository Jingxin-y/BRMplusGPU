data.generation <- function(param, n, alpha.true, beta.true, gamma.true){
  
  getProb = if (param == "RR") getProbRR else getProbRD
  
  v.1         = rep(1,n)       # intercept term
  v.2         = runif(n,0,0.6)
  v           = cbind(v.1,v.2)
  v.1 = as.matrix(v.1, ncol = 1)
  pscore.true = exp(v %*% gamma.true) / (1+exp(v %*% gamma.true))
  p0p1.true   = getProb(v %*% alpha.true,v %*% beta.true)
  x           = rbinom(n, 1, pscore.true)
  pA.true       = p0p1.true[,1]
  pA.true[x==1] = p0p1.true[x==1,2]
  y = rbinom(n, 1, pA.true)
  
  summary(predict(glm(x~v.1,family=binomial(link = "log")),type = "response"))
  
  Na0 <- sum(x==0)
  Na1 <- sum(x==1)
  N0_1 <- sum(y[which(x==0)])
  N1_1 <- sum(y[which(x==1)])
  
  data.simulation <- list(data = data.frame(y,x,v), count = c(Na0,Na1,N0_1,N1_1))
  return(data.simulation)
}

simulate.rr <- function(n, event, hypothesis){
  
  if (event == "common"){
    if (hypothesis == "null"){
      alpha.true <- c(0,0)
      beta.true  <- c(1.5, 0.6)
      gamma.true <- c(0.2, -0.5)
    }else{
      alpha.true <- c(0.15,0.5)
      beta.true  <- c(1.65, 0.5)
      gamma.true <- c(0.2, -0.5)
    }
  }else{
    if (hypothesis == "null"){
      alpha.true <- c(0,0)
      beta.true  <- c(-4.7, 0.5)
      gamma.true <- c(0.2, -0.5)
    }else{
      alpha.true <- c(0.55,0.5)
      beta.true  <- c(-5.5, 0.5)
      gamma.true <- c(0.2, -0.5)
    }
  }
  
  data.simulation <- data.generation('RR', n, alpha.true, beta.true, gamma.true)
  
  va = cbind(data.simulation$data$v.1,data.simulation$data$v.2)
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
  
  
  
  ##brm+exact
  est.exact <- exact('RR', y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-3, est.brm$point.est, est.brm$se.est, pa, pb)
  # est.exact.ad <- exact('RR', y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-3, est.brm.ad$point.est, est.brm.ad$se.est, pa, pb)
  
  ###result
  point.est <- c(est.brm$point.est[1:2],est.brm$point.est[1:2])
  se.est <- c(est.brm$se.est[1:2],est.brm$se.est[1:2])
  con.lower <- c(est.brm$conf.lower[1:2],est.exact$low)
  con.upper <- c(est.brm$conf.upper[1:2],est.exact$up)
  p.value <- c(est.brm$p.value[1:2],est.exact$p)
  
  result.comp <- rbind(point.est,se.est,con.lower,con.upper,p.value)
  colnames(result.comp) <- c("brm_alpha0","brm_alpha1",
                             "brm_exact_alpha0","brm_exact_alpha1")
  return(result.comp)
}

run <- function(param,n,event,hypothesis){
  simulate.fun = if (param == "RR") simulate.rr else simulate.rd
  result = simulate.fun(n,event,hypothesis)
  return(result)
}