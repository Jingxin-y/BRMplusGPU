#' Generate Simulated Binary Data under RR/RD Models
#'
#' @description
#' Simulates a binary outcome \code{y} with a binary treatment \code{x} and two
#' covariates \code{v.1} (intercept) and \code{v.2} (uniform on \eqn{[0,0.6]}).
#' Treatment assignment follows a logistic propensity score with linear predictor
#' \code{v %*% gamma.true}. Outcome probabilities under \code{x=0,1} are produced
#' by \code{getProbRR} (for \code{param="RR"}) or \code{getProbRD} (for
#' \code{param="RD"}), which must return a two-column matrix \code{cbind(p0, p1)}.
#'
#' @param param Character. \code{"RR"} (relative risk) or \code{"RD"} (risk difference).
#' @param n Integer. Sample size.
#' @param alpha.true Numeric vector (\eqn{p_a \times 1}). Structural parameter(s) for \eqn{\alpha}.
#' @param beta.true Numeric vector (\eqn{p_b \times 1}). Structural parameter(s) for \eqn{\beta}.
#' @param gamma.true Numeric vector (length 2). Coefficients for treatment model (\code{pscore}).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{data}}{A \code{data.frame} with columns \code{y}, \code{x}, \code{v.1}, \code{v.2}.}
#'   \item{\code{count}}{Numeric vector \code{c(Na0, Na1, N0_1, N1_1)} giving group sizes
#'   and number of successes by arm.}
#' }
#'

data.generation <- function(param, n, alpha.true, beta.true, gamma.true){

  getProb = if (param == "RR") getProbRR else getProbRD

  v.1         = rep(1,n)       # intercept term
  v.2         = runif(n,0,0.6)
  v           = cbind(v.1,v.2)
  v.1 = as.matrix(v.1, ncol = 1)
  pscore.true = exp(v %*% gamma.true) / (1+exp(v %*% gamma.true))
  p0p1.true   = getProb(v.1 %*% alpha.true,v %*% beta.true)
  x           = rbinom(n, 1, pscore.true)
  pA.true       = p0p1.true[,1]
  pA.true[x==1] = p0p1.true[x==1,2]
  y = rbinom(n, 1, pA.true)

  Na0 <- sum(x==0)
  Na1 <- sum(x==1)
  N0_1 <- sum(y[which(x==0)])
  N1_1 <- sum(y[which(x==1)])

  data.simulation <- list(data = data.frame(y,x,v), count = c(Na0,Na1,N0_1,N1_1))
  return(data.simulation)
}


#' Quasi-Poisson Log-Link with Robust (HC0) SE for Treatment Effect
#'
#' @description
#' Fits \code{glm(y ~ x + v.1 + v.2 - 1, family = quasipoisson(link="log"))} and
#' reports the coefficient, robust standard error (HC0), Wald CI, and two-sided
#' p-value for the \code{x} effect (assumed to be the first coefficient).
#'
#' @param data \code{data.frame} with columns \code{y}, \code{x}, \code{v.1}, \code{v.2}.
#'
#' @return Numeric vector \code{c(est, se.robust, lower, upper, p.value)} for the
#' treatment coefficient on the log scale.
#'
#' @details
#' The robust variance uses \code{vcovHC(fit, type="HC0")}. This can be viewed as
#' a log–Poisson working model with overdispersion and sandwich SEs.


quasi.poisson <- function(data){
  fit.qp <- glm(y~x+v.1+v.2-1, family = quasipoisson(link = "log"), data = data)
  vc <- vcovHC(fit.qp, type = "HC0")[1,1]
  est <- coef(fit.qp)[1]
  se.robust  <- sqrt(vc)
  p.robust   <- 2 * min(pnorm(abs(est/se.robust)),(1 - pnorm(abs(est/se.robust))))
  lower <- est - 1.96 * se.robust
  upper <- est + 1.96 * se.robust

  return(c(est,se.robust,lower,upper,p.robust))
}

#' Simulate and Compare RR Estimators Across Multiple Methods
#'
#' @description
#' Generates data under an RR parametrization and computes estimates, SEs, CIs,
#' and p-values for a suite of methods: BRM MLE, BRM+adaptive (Bayes fallback),
#' CMH, log-binomial, log-Poisson, robust log-Poisson (quasi-Poisson + HC0),
#' BRM+Firth, profile-exact (based on BRM), and several g-computation variants
#' (plain, bias-reduction, Firth-corrected/FC with bias-reduction BR1/BR2). Returns a 5×14 matrix:
#' rows = \code{point.est}, \code{se.est}, \code{con.lower}, \code{con.upper}, \code{p.value};
#' columns labeled by method.
#'
#' @param n Integer. Sample size.
#' @param event Character. \code{"common"} or \code{"rare"} to set truth.
#' @param hypothesis Character. \code{"null"} or \code{"alternative"}.
#'
#' @return A numeric matrix with rows \code{point.est}, \code{se.est},
#' \code{con.lower}, \code{con.upper}, \code{p.value} and 14 method columns:
#' \code{c("brm","brm_ad","CMH","log-binomial","log-poisson","robust log-possion",
#' "brm_firth","brm_exact","brm_exact_ad","g-computation","GC_BR","GC_FC","GC_FC_BR1","GC_FC_BR2")}.
#'

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
  est.lb <- glm(y~x+v.1+v.2-1, family = binomial(link = "log"), data = data.simulation$data, start = rep(-0.01,3))

  ##log-poisson
  est.lp <- glm(y~x+v.1+v.2-1, family = poisson(link = "log"), data = data.simulation$data)

  ##robust log-poisson

  est.rlp <- quasi.poisson(data.simulation$data)

  ##brm + firth

  est.brm.firth <- MLEst('RR', y, x, va, vb, weight, max.step, thres, alpha.start = rep(0, pa),
                               beta.start = rep(0, pb), pa, pb, method="firth")

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

  ##g-computaion & g-computation_BR
  Y1 <- y[which(x==1)]
  Y0 <- y[which(x==0)]
  V2.1 <- v.2[which(x==1)]
  V2.0 <- v.2[which(x==0)]
  X1 <- x[which(x==1)]
  X0 <- x[which(x==0)]

  data.treat <- data.frame(Y1,V2.1)
  data.control <- data.frame(Y0,V2.0)

  est.treat <- glm(Y1~V2.1, family = binomial, data = data.treat)
  est.control <- glm(Y0~V2.0, family = binomial, data = data.control)

  beta.hat.treat <- est.treat$coefficients
  beta.hat.control <- est.control$coefficients

  V.FC.treat <- cbind(1,V2.1)
  V.FC.control <- cbind(1,V2.0)


  beta.hat.star.treat <- beta.hat.treat + colMeans(hatvalues(est.treat)*phi(Y1,V.FC.treat,beta.hat.treat,sum(x==1)/n))
  beta.hat.star.control <- beta.hat.control + colMeans(hatvalues(est.control)*phi(Y0,V.FC.control,beta.hat.control,sum(x==0)/n))

  #beta_hat
  p.hat.treat <- mean(c(Y1,predict(est.treat,newdata = data.control, type = "response")))
  p.hat.control <- mean(c(Y0,predict(est.control,newdata = data.treat, type = "response")))
  alpha.hat <- log(p.hat.treat/p.hat.control)

  #beta_hat_star
  p.hat.star.treat <-mean(c(Y1,m(V.FC.control%*%beta.hat.star.treat)))
  p.hat.star.control <- mean(c(Y0,m(V.FC.treat%*%beta.hat.star.control)))
  alpha.hat.star <- log(p.hat.star.treat/p.hat.star.control)

  li.hat <- l.mu(Y1,V.FC.treat,beta.hat.treat,Y0,V.FC.control,beta.hat.control)
  li.hat.star <- l.mu(Y1,V.FC.treat,beta.hat.star.treat,Y0,V.FC.control,beta.hat.star.control)

  se.hat <- sqrt(var.est.RR(li.hat,p.hat.control,p.hat.treat))
  se.hat.star <- sqrt(var.est.RR(li.hat.star,p.hat.star.control,p.hat.star.treat))

  #beta_tilde
  fit.treat <- logistf(Y1 ~ V2.1,data = data.treat)
  fit.control <- logistf(Y0 ~ V2.0,data = data.control)

  beta.tilde.treat <- fit.treat$coefficients
  beta.tilde.control <- fit.control$coefficients

  #beta_tilde_star
  beta.tilde.star.treat <- beta.tilde.treat + colMeans(as.vector(hii(V.FC.treat,beta.tilde.treat))*(phi(Y1,V.FC.treat,beta.tilde.treat,sum(x==1)/n)
                                                                                                    -(V.FC.treat*as.vector(1-2*m(V.FC.treat%*%beta.tilde.treat)))%*%t(ginv(fish(V.FC.treat,beta.tilde.treat)))/2))
  beta.tilde.star.control <- beta.tilde.control + colMeans(as.vector(hii(V.FC.control,beta.tilde.control))*(phi(Y0,V.FC.control,beta.tilde.control,sum(x==0)/n)
                                                                                                            -(V.FC.control*as.vector(1-2*m(V.FC.control%*%beta.tilde.control)))%*%t(ginv(fish(V.FC.control,beta.tilde.control)))/2))
  #beta_tilde_doustar
  beta.tilde.doustar.treat <- beta.tilde.treat - colMeans(as.vector(hii(V.FC.treat,beta.tilde.treat))*((V.FC.treat*as.vector(1-2*m(V.FC.treat%*%beta.tilde.treat)))%*%t(ginv(fish(V.FC.treat,beta.tilde.treat)))/2))
  beta.tilde.doustar.control <- beta.tilde.control - colMeans(as.vector(hii(V.FC.control,beta.tilde.control))*((V.FC.control*as.vector(1-2*m(V.FC.control%*%beta.tilde.control)))%*%t(ginv(fish(V.FC.control,beta.tilde.control)))/2))

  #beta_tilde
  p.tilde.treat <-mean(c(Y1,m(V.FC.control%*%beta.tilde.treat)))
  p.tilde.control <- mean(c(Y0,m(V.FC.treat%*%beta.tilde.control)))
  alpha.tilde <- log(p.tilde.treat/p.tilde.control)

  #beta_tilde_star
  p.tilde.star.treat <-mean(c(Y1,m(V.FC.control%*%beta.tilde.star.treat)))
  p.tilde.star.control <- mean(c(Y0,m(V.FC.treat%*%beta.tilde.star.control)))
  alpha.tilde.star <- log(p.tilde.star.treat/p.tilde.star.control)

  #beta_tilde_starstar
  p.tilde.doustar.treat <-mean(c(Y1,m(V.FC.control%*%beta.tilde.doustar.treat)))
  p.tilde.doustar.control <- mean(c(Y0,m(V.FC.treat%*%beta.tilde.doustar.control)))
  alpha.tilde.doustar <- log(p.tilde.doustar.treat/p.tilde.doustar.control)

  li.tilde <- l.mu(Y1,V.FC.treat,beta.tilde.treat,Y0,V.FC.control,beta.tilde.control)
  li.tilde.star <- l.mu(Y1,V.FC.treat,beta.tilde.star.treat,Y0,V.FC.control,beta.tilde.star.control)
  li.tilde.doustar <- l.mu(Y1,V.FC.treat,beta.tilde.doustar.treat,Y0,V.FC.control,beta.tilde.doustar.control)

  se.tilde <- sqrt(var.est.RR(li.tilde,p.tilde.control,p.tilde.treat))
  se.tilde.star <- sqrt(var.est.RR(li.tilde.star,p.tilde.star.control,p.tilde.star.treat))
  se.tilde.doustar <- sqrt(var.est.RR(li.tilde.doustar,p.tilde.doustar.control,p.tilde.doustar.treat))



  ##brm+exact
  est.exact <- exact('RR', y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-3, est.brm$point.est, est.brm$se.est, pa, pb)
  est.exact.ad <- exact('RR', y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-3, est.brm.ad$point.est, est.brm.ad$se.est, pa, pb)

  ###result
  point.est <- as.vector(c(est.brm$point.est[1],
                 est.brm.ad$point.est[1],
                 log(est.CMH$measure[2,1]),
                 est.lb$coefficients[1],
                 est.lp$coefficients[1],
                 est.rlp[1],
                 est.brm.firth$point.est[1],
                 est.brm$point.est[1],
                 est.brm.ad$point.est[1],
                 alpha.hat,
                 alpha.hat.star,
                 alpha.tilde,
                 alpha.tilde.star,
                 alpha.tilde.doustar))
  se.est <- as.vector(c(est.brm$se.est[1],
              est.brm.ad$se.est[1],
              (log(est.CMH$measure[2,1])-log(est.CMH$measure[2,2]))/qnorm(0.975),
              summary(est.lb)$coefficients[1,2],
              summary(est.lp)$coefficients[1,2],
              est.rlp[2],
              est.brm.firth$se.est[1],
              est.brm$se.est[1],
              est.brm.ad$se.est[1],
              se.hat,
              se.hat.star,
              se.tilde,
              se.tilde.star,
              se.tilde.doustar))
  con.lower <- as.vector(c(est.brm$conf.lower[1],
                 est.brm.ad$conf.lower[1],
                 log(est.CMH$measure[2,2]),
                 confint.default(est.lb,level = 0.95)[1,1],
                 confint.default(est.lp,level = 0.95)[1,1],
                 est.rlp[3],
                 est.brm.firth$conf.lower[1],
                 est.exact$low[1],
                 est.exact.ad$low[1],
                 alpha.hat-qnorm(0.975)*se.hat,
                 alpha.hat.star-qnorm(0.975)*se.hat.star,
                 alpha.tilde-qnorm(0.975)*se.tilde,
                 alpha.tilde.star-qnorm(0.975)*se.tilde.star,
                 alpha.tilde.doustar-qnorm(0.975)*se.tilde.doustar))
  con.upper <- as.vector(c(est.brm$conf.upper[1],
                 est.brm.ad$conf.upper[1],
                 log(est.CMH$measure[2,3]),
                 confint.default(est.lb,level = 0.95)[1,2],
                 confint.default(est.lp,level = 0.95)[1,2],
                 est.rlp[4],
                 est.brm.firth$conf.upper[1],
                 est.exact$up[1],
                 est.exact.ad$up[1],
                 alpha.hat+qnorm(0.975)*se.hat,
                 alpha.hat.star+qnorm(0.975)*se.hat.star,
                 alpha.tilde+qnorm(0.975)*se.tilde,
                 alpha.tilde.star+qnorm(0.975)*se.tilde.star,
                 alpha.tilde.doustar+qnorm(0.975)*se.tilde.doustar))

  p.value <- as.vector(c(est.brm$p.value[1],
               est.brm.ad$p.value[1],
               est.CMH$p.value[2,1],
               summary(est.lb)$coefficients[1,4],
               summary(est.lp)$coefficients[1,4],
               est.rlp[5],
               est.brm.firth$p.value[1],
               est.exact$p[1],
               est.exact.ad$p[1],
               2*min(pnorm(alpha.hat/se.hat),1-pnorm(alpha.hat/se.hat)),
               2*min(pnorm(alpha.hat.star/se.hat.star),1-pnorm(alpha.hat.star/se.hat.star)),
               2*min(pnorm(alpha.tilde/se.tilde),1-pnorm(alpha.tilde/se.tilde)),
               2*min(pnorm(alpha.tilde.star/se.tilde.star),1-pnorm(alpha.tilde.star/se.tilde.star)),
               2*min(pnorm(alpha.tilde.doustar/se.tilde.doustar),1-pnorm(alpha.tilde.doustar/se.tilde.doustar))))

  result.comp <- rbind(point.est,se.est,con.lower,con.upper,p.value)
  colnames(result.comp) <- c("brm","brm_ad","CMH","log-binomial","log-poisson","robust log-possion","brm_firth",
                             "brm_exact","brm_exact_ad","g-computation","GC_BR","GC_FC","GC_FC_BR1","GC_FC_BR2")
  return(result.comp)
}

#' Simulate and Compare RD Estimators Across Multiple Methods
#'
#' @description
#' Generates data under an RD parametrization and computes estimates, SEs, CIs,
#' and p-values for multiple methods: BRM MLE (original and adaptive Bayes),
#' Bayesian RD with simple conjugate prior, GLM with identity link (if feasible),
#' LPM with robust SEs, Miettinen–Nurminen (MN), BRM+Firth, profile-exact, and
#' g-computation variants (plain, BR, FC and BR1/BR2). Returns a 5×14 matrix with
#' rows \code{point.est}, \code{se.est}, \code{CI.low.or}, \code{CI.up.or}, \code{p.value}.
#'
#' @param n Integer. Sample size.
#' @param event Character. \code{"common"} or \code{"rare"} to set truth.
#' @param hypothesis Character. \code{"null"} or \code{"alternative"}.
#'
#' @return A numeric matrix with rows \code{point.est}, \code{se.est},
#' \code{CI.low.or}, \code{CI.up.or}, \code{p.value} and 14 method columns:
#' \code{c("brm","brm_ad","bayes","glm","lpm","MN","firth","brm_exact","brm_exact_ad",
#' "g-computation","GC_BR","GC_FC","GC_FC_BR1","GC_FC_BR2")}.

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
  e.glm <- glm(y~x+v.1+v.2-1, family = binomial(link = "identity"), data = data.simulation$data,start = rep(0.01,3))
  est.glm <- get_estimate(e.glm$coefficients[1], summary(e.glm)$coefficients[1,2], as.numeric(confint.default(e.glm,level = 0.95)[1,]))

  ## Linear probability model (LPM) + robust SE
  lpm <- lm(y~x+v.1+v.2-1,data = data.simulation$data)
  e.lpm <- coeftest(lpm, vcov = vcovHC(lpm,type = "HC3"))
  est.lpm <- get_estimate(e.lpm[1,1],e.lpm[1,2],as.numeric(c(e.lpm[1,1]-1.96*e.lpm[1,2],e.lpm[1,1]+1.96*e.lpm[1,2])))

  ## Miettinen–Nurminen
  est.MN.point <- P1-P0
  est.MN.CI <- diffscoreci(N1_1, Na1, N0_1, Na0, conf.level = 0.95)
  est.MN.se <- (est.MN.CI$conf.int[2]-est.MN.CI$conf.int[1])/(2*qnorm(0.975))
  est.MN <- get_estimate(est.MN.point,est.MN.se,c(est.MN.CI$conf.int[1],est.MN.CI$conf.int[2]))
  #  p.MN <- 2*(1-pnorm(abs(z2stat(N1_1,Na1,N0_1,Na0,dif=0))))


  ##g-computaion & g-computation_BR
  Y1 <- y[which(x==1)]
  Y0 <- y[which(x==0)]
  V2.1 <- v.2[which(x==1)]
  V2.0 <- v.2[which(x==0)]
  X1 <- x[which(x==1)]
  X0 <- x[which(x==0)]

  data.treat <- data.frame(Y1,V2.1)
  data.control <- data.frame(Y0,V2.0)

  est.treat <- glm(Y1~V2.1, family = binomial, data = data.treat)
  est.control <- glm(Y0~V2.0, family = binomial, data = data.control)

  beta.hat.treat <- est.treat$coefficients
  beta.hat.control <- est.control$coefficients

  V.FC.treat <- cbind(1,V2.1)
  V.FC.control <- cbind(1,V2.0)


  beta.hat.star.treat <- beta.hat.treat + colMeans(hatvalues(est.treat)*phi(Y1,V.FC.treat,beta.hat.treat,sum(x==1)/n))
  beta.hat.star.control <- beta.hat.control + colMeans(hatvalues(est.control)*phi(Y0,V.FC.control,beta.hat.control,sum(x==0)/n))

  #beta_tilde
  fit.treat <- logistf(Y1 ~ V2.1,data = data.treat)
  fit.control <- logistf(Y0 ~ V2.0,data = data.control)

  beta.tilde.treat <- fit.treat$coefficients
  beta.tilde.control <- fit.control$coefficients

  #beta_tilde_star
  beta.tilde.star.treat <- beta.tilde.treat + colMeans(as.vector(hii(V.FC.treat,beta.tilde.treat))*(phi(Y1,V.FC.treat,beta.tilde.treat,sum(x==1)/n)
                                                                                                    -(V.FC.treat*as.vector(1-2*m(V.FC.treat%*%beta.tilde.treat)))%*%t(ginv(fish(V.FC.treat,beta.tilde.treat)))/2))
  beta.tilde.star.control <- beta.tilde.control + colMeans(as.vector(hii(V.FC.control,beta.tilde.control))*(phi(Y0,V.FC.control,beta.tilde.control,sum(x==0)/n)
                                                                                                            -(V.FC.control*as.vector(1-2*m(V.FC.control%*%beta.tilde.control)))%*%t(ginv(fish(V.FC.control,beta.tilde.control)))/2))
  #beta_tilde_doustar
  beta.tilde.doustar.treat <- beta.tilde.treat - colMeans(as.vector(hii(V.FC.treat,beta.tilde.treat))*((V.FC.treat*as.vector(1-2*m(V.FC.treat%*%beta.tilde.treat)))%*%t(ginv(fish(V.FC.treat,beta.tilde.treat)))/2))
  beta.tilde.doustar.control <- beta.tilde.control - colMeans(as.vector(hii(V.FC.control,beta.tilde.control))*((V.FC.control*as.vector(1-2*m(V.FC.control%*%beta.tilde.control)))%*%t(ginv(fish(V.FC.control,beta.tilde.control)))/2))

  p.hat.treat <- mean(c(Y1,m(V.FC.control%*%beta.hat.treat)))
  p.hat.control <- mean(c(Y0,m(V.FC.treat%*%beta.hat.control)))
  alpha.hat <- atanh(p.hat.treat-p.hat.control)

  #beta_hat_star
  p.hat.star.treat <-mean(c(Y1,m(V.FC.control%*%beta.hat.star.treat)))
  p.hat.star.control <- mean(c(Y0,m(V.FC.treat%*%beta.hat.star.control)))
  alpha.hat.star <- atanh(p.hat.star.treat-p.hat.star.control)

  #beta_tilde
  p.tilde.treat <-mean(c(Y1,m(V.FC.control%*%beta.tilde.treat)))
  p.tilde.control <- mean(c(Y0,m(V.FC.treat%*%beta.tilde.control)))
  alpha.tilde <- atanh(p.tilde.treat-p.tilde.control)

  p.tilde.star.treat <-mean(c(Y1,m(V.FC.control%*%beta.tilde.star.treat)))
  p.tilde.star.control <- mean(c(Y0,m(V.FC.treat%*%beta.tilde.star.control)))
  alpha.tilde.star <- atanh(p.tilde.star.treat-p.tilde.star.control)

  #beta_tilde_starstar
  p.tilde.doustar.treat <-mean(c(Y1,m(V.FC.control%*%beta.tilde.doustar.treat)))
  p.tilde.doustar.control <- mean(c(Y0,m(V.FC.treat%*%beta.tilde.doustar.control)))
  alpha.tilde.doustar <- atanh(p.tilde.doustar.treat-p.tilde.doustar.control)


  li.hat <- l.mu(Y1,V.FC.treat,beta.hat.treat,Y0,V.FC.control,beta.hat.control)
  li.hat.star <- l.mu(Y1,V.FC.treat,beta.hat.star.treat,Y0,V.FC.control,beta.hat.star.control)
  li.tilde <- l.mu(Y1,V.FC.treat,beta.tilde.treat,Y0,V.FC.control,beta.tilde.control)
  # li.firth <- l.mu(Y1,V.FC.treat,beta.firth.treat,Y0,V.FC.control,beta.firth.control)
  li.tilde.star <- l.mu(Y1,V.FC.treat,beta.tilde.star.treat,Y0,V.FC.control,beta.tilde.star.control)
  li.tilde.doustar <- l.mu(Y1,V.FC.treat,beta.tilde.doustar.treat,Y0,V.FC.control,beta.tilde.doustar.control)

  se.hat <- sqrt(var.est.RD(li.hat,p.hat.control,p.hat.treat))
  se.hat.star <- sqrt(var.est.RD(li.hat.star,p.hat.star.control,p.hat.star.treat))
  se.tilde <- sqrt(var.est.RD(li.tilde,p.tilde.control,p.tilde.treat))
  # se.firth <- sqrt(var.est(li.firth,p.firth.control,p.firth.treat))
  se.tilde.star <- sqrt(var.est.RD(li.tilde.star,p.tilde.star.control,p.tilde.star.treat))
  se.tilde.doustar <- sqrt(var.est.RD(li.tilde.doustar,p.tilde.doustar.control,p.tilde.doustar.treat))


  ##brm_firth
  est.brm.Firth <- MLEst('RD', y, x, va, vb, weight, max.step, thres, alpha.start = rep(0, pa),
  beta.start = rep(0, pb), pa, pb, method="firth")



  #### CI and p.value

  est.exact.ad <- exact('RD', y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-3, est.brm.ad$point.est, est.brm.ad$se.est, pa, pb)
  est.exact <- exact('RD', y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-3, est.brm.or$point.est, est.brm.or$se.est, pa, pb)


  ###result
  point.est <- c(est.brm.or$point.est[1],
                 est.brm.ad$point.est[1],
                 est.bayes$point.est,
                 est.glm$point.est,
                 est.lpm$point.est,
                 est.MN$point.est,
                 est.brm.Firth$point.est[1],
                 est.brm.or$point.est[1],
                 est.brm.ad$point.est[1],
                 alpha.hat,
                 alpha.hat.star,
                 alpha.tilde,
                 alpha.tilde.star,
                 alpha.tilde.doustar)
  se.est <- c(est.brm.or$se.est[1],
              est.brm.ad$se.est[1],
              est.bayes$se.est,
              est.glm$se.est,
              est.lpm$se.est,
              est.MN$se.est,
              est.brm.Firth$se.est[1],
              est.brm.or$se.est[1],
              est.brm.ad$se.est[1],
              se.hat,
              se.hat.star,
              se.tilde,
              se.tilde.star,
              se.tilde.doustar)
  CI.low.or <- c(est.brm.or$conf.lower[1],
                 est.brm.ad$conf.lower[1],
                 est.bayes$conf.lower,
                 est.glm$CI[1],
                 est.lpm$CI[1],
                 est.MN$CI[1],
                 est.brm.Firth$conf.lower[1],
                 est.exact$low[1],
                 est.exact.ad$low[1],
                 alpha.hat-qnorm(0.975)*se.hat,
                 alpha.hat.star-qnorm(0.975)*se.hat.star,
                 alpha.tilde-qnorm(0.975)*se.tilde,
                 alpha.tilde.star-qnorm(0.975)*se.tilde.star,
                 alpha.tilde.doustar-qnorm(0.975)*se.tilde.doustar)
  CI.up.or <- c(est.brm.or$conf.upper[1],
                est.brm.ad$conf.upper[1],
                est.bayes$conf.upper,
                est.glm$CI[2],
                est.lpm$CI[2],
                est.MN$CI[2],
                est.brm.Firth$conf.upper[1],
                est.exact$up[1],
                est.exact.ad$up[1],
                alpha.hat+qnorm(0.975)*se.hat,
                alpha.hat.star+qnorm(0.975)*se.hat.star,
                alpha.tilde+qnorm(0.975)*se.tilde,
                alpha.tilde.star+qnorm(0.975)*se.tilde.star,
                alpha.tilde.doustar+qnorm(0.975)*se.tilde.doustar)
  p.value <- c(est.brm.or$p.value[1],
               est.brm.ad$p.value[1],
               est.bayes$p.value,
               summary(e.glm)$coefficients[1,4],
               summary(lpm)$coefficients[1,4],
               min(pnorm(alpha.hat/se.hat),1-pnorm(alpha.hat/se.hat)),
               est.brm.Firth$p.value[1],
               est.exact$p[1],
               est.exact.ad$p[1],
               2*min(pnorm(est.MN$point.est/est.MN$se.est),1-pnorm(est.MN$point.est/est.MN$se.est)),
               2*min(pnorm(alpha.hat.star/se.hat.star),1-pnorm(alpha.hat.star/se.hat.star)),
               2*min(pnorm(alpha.tilde/se.tilde),1-pnorm(alpha.tilde/se.tilde)),
               2*min(pnorm(alpha.tilde.star/se.tilde.star),1-pnorm(alpha.tilde.star/se.tilde.star)),
               2*min(pnorm(alpha.tilde.doustar/se.tilde.doustar),1-pnorm(alpha.tilde.doustar/se.tilde.doustar)))

  result.comp <- rbind(point.est,se.est,CI.low.or,CI.up.or,p.value)
  colnames(result.comp) <- c("brm","brm_ad","bayes","glm","lpm","MN", "firth",
                             "brm_exact","brm_exact_ad","g-computation","GC_BR","GC_FC","GC_FC_BR1","GC_FC_BR2")
  return(result.comp)
}

#' Run a Single Simulation for RR or RD
#'
#' @description
#' Dispatch helper that runs \code{simulate.rr()} if \code{param="RR"} and
#' \code{simulate.rd()} otherwise.
#'
#' @param param Character. \code{"RR"} or \code{"RD"}.
#' @param n Integer. Sample size.
#' @param event Character. \code{"common"} or \code{"rare"}.
#' @param hypothesis Character. \code{"null"} or \code{"alternative"}.
#'
#' @return The matrix returned by the corresponding simulator.


run <- function(param,n,event,hypothesis){
  simulate.fun = if (param == "RR") simulate.rr else simulate.rd
  result = simulate.fun(n,event,hypothesis)
  return(result)
}

