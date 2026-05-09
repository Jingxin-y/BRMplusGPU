#' Maximum Likelihood Estimation for Risk‐Difference models with Firth's Augmentation
#'
#' Firth's method: Firth, D. (1993). Bias reduction of maximum likelihood estimates. Biometrika, 80(1), 27-38.
#'
#' Optimize the log‐likelihood for a binary‐outcome regression model on the
#' risk‐difference (RD) scale using an iterative augmentation scheme.
#'
#' @param param Character string, takes\code{"RD"}
#'
#' @param y Numeric vector of length \eqn{n}.  Binary outcome values (0/1).
#' @param x Numeric vector of length \eqn{n}.  Binary exposure indicator (0/1).
#' @param va Numeric matrix of dimension \eqn{n \times pa}.
#' @param vb Numeric matrix of dimension \eqn{n \times pb}.
#' @param alpha.start Numeric vector of length \eqn{pa}, or \code{NULL}.  Initial
#'   values for the \eqn{\alpha} parameters; if \code{NULL}, defaults to a zero
#'   vector of length \eqn{pa}.
#' @param beta.start Numeric vector of length \eqn{pb}, or \code{NULL}.  Initial
#'   values for the \eqn{\beta} parameters; if \code{NULL}, defaults to a zero
#'   vector of length \eqn{pb}.
#' @param weight Numeric vector of length \eqn{n}.
#' @param max.step Integer.  Maximum number of alternating iterations to perform.
#' @param thres Numeric.  Convergence threshold on relative change in parameters.
#' @param pa Integer.  Number of \eqn{\alpha} parameters (\eqn{pa}).
#' @param pb Integer.  Number of \eqn{\beta} parameters (\eqn{pb}).
#'
#' this model is for pa = 1 and pb = 2, if pa and pb change, or param change, the function
#' 'compute.components', 'compute.augmentation' and 'compute.score' need to be changed.
#' If model change, 'compute.components' and 'compute.score' need to be changed.
#'

max.likelihood.firth.rd = function(param, y, x, va, vb, alpha.start, beta.start, weight, max.step, thres, pa, pb) {
  ### augmentation calculation, compute:
  #'   - the observed Fisher information matrix,
  #'   - its inverse,
  #'   - third‐order cumulants \(k_{s t u}\),  κ_{s,t,u} = n^{-1} * E{ U_s, U_t, U_u }
  #'   - mixed cumulants \(k_{s, t u}\), κ_{s,tu}  = n^{-1} * E{ U_s, U_{tu}}.
  #'
compute.components = function(x, alpha.ml, beta.ml, va, vb, weight) {


  p0p1 = getProbRD(va %*% alpha.ml, vb %*% beta.ml) #n by 2
  # p0p1 = cbind(p0, p1): n * 2 matrix
  p0 = p0p1[, 1]
  p1 = p0p1[, 2]
  n = nrow(vb)
  pA = p0
  pA[x == 1] = p1[x == 1]
  # s0 = p0 * (1 - p0) # n by 1
  # s1 = p1 * (1 - p1)
  # sA = pA * (1 - pA)
  s0 <- pmax(p0 * (1 - p0), 1e-8)
  s1 <- pmax(p1 * (1 - p1), 1e-8)
  sA <- pmax(pA * (1 - pA), 1e-8)


  rho = as.vector(tanh(va %*% alpha.ml))  #estimated risk differences n by 1

  ### First order derivatives ###
  va.1 = 1
  vb.1 = 1
  expect.dl.by.dpA = 1/sA # n by 1
  dp0.by.dphi = s0 * s1/(s0 + s1) # n by 1
  dp0.by.drho = -s0/(s0 + s1) # n by 1
  drho.by.dalpha = va.1*(1 - rho^2) # n by 1
  dphi.by.dbeta =vb.1

  dpA.by.drho = dp0.by.drho + x # n by 1
  dpA.by.dalpha = drho.by.dalpha * dpA.by.drho # n by 1
  dpA.by.dphi = dp0.by.dphi # n by 1
  dpA.by.dbeta = dphi.by.dbeta * dpA.by.dphi # n by 1




  ### Second order derivatives ###

  expect.d2l.by.dpA.2 = -(1 - 2*pA)/sA^2 # n by 1
  d2pA.by.drho.2 = s0 * s1 * (2 - 2 * p0 - 2 * p1)/(s0 + s1)^3 # n by 1
  d2pA.by.dphi.drho = (s0 * (1 - 2 * p1) - s1 * (1 - 2 * p0)) * s0 * s1/(s0 + s1)^3 # n by 1
  d2pA.by.dphi.2 = (s0^2 * (1 - 2 * p1) + s1^2 * (1 - 2 * p0)) * s0 * s1/(s0 +s1)^3 # n by 1

  d2rho.by.dalpha.2 = -2 * va.1 * rho * drho.by.dalpha # n by 1

  ### Compute elements of the Hessian matrix ###

  d2l.by.dalpha.2 = dpA.by.dalpha * expect.d2l.by.dpA.2 * weight * dpA.by.dalpha +
    drho.by.dalpha * expect.dl.by.dpA * d2pA.by.drho.2 * weight * drho.by.dalpha -
    2 * va.1 * rho * expect.dl.by.dpA * dpA.by.drho * weight * drho.by.dalpha # n by 1

  d2l.by.dalpha.dbeta = dpA.by.dalpha * expect.d2l.by.dpA.2 * weight * dpA.by.dbeta +
    drho.by.dalpha * expect.dl.by.dpA * d2pA.by.dphi.drho * weight * dphi.by.dbeta
  d2l.by.dbeta.dalpha = d2l.by.dalpha.dbeta

  d2l.by.dbeta.2 = dpA.by.dbeta * expect.d2l.by.dpA.2 * weight * dpA.by.dbeta +
    dphi.by.dbeta * expect.dl.by.dpA * d2pA.by.dphi.2 * weight * dphi.by.dbeta

  ###


  ##  fisher info

  ## k_{s,t,u}
  c.stu.A = (1-2*pA)/(sA^2)
  c.stu.alpha = (dp0.by.drho + x) * drho.by.dalpha
  c.stu.beta = dp0.by.dphi

  k.aaa = c.stu.A*c.stu.alpha^3
  k.aab = c.stu.A*c.stu.alpha^2*c.stu.beta
  k.abb = c.stu.A*c.stu.alpha*c.stu.beta^2
  k.bbb = c.stu.A*c.stu.beta^3


  ## k_{s,tu}


  k.a.aa = c.stu.alpha*d2l.by.dalpha.2
  k.a.ab = c.stu.alpha*d2l.by.dalpha.dbeta
  k.a.bb = c.stu.alpha*d2l.by.dbeta.2
  k.b.aa = c.stu.beta*d2l.by.dalpha.2
  k.b.ab = c.stu.beta*d2l.by.dalpha.dbeta
  k.b.bb = c.stu.beta*d2l.by.dbeta.2

  expect.dl.by.dpA.squared = 1/sA
  dp0.by.dphi = s0 * s1/(s0 + s1)
  dp0.by.drho = -s0/(s0 + s1)
  drho.by.dalpha = va*(1 - rho^2)
  dphi.by.dbeta = vb

  tmp = cbind((dp0.by.drho + x) * drho.by.dalpha, dp0.by.dphi * dphi.by.dbeta)
  fisher.info = (t(expect.dl.by.dpA.squared * weight * tmp) %*% tmp)


  return(list(fisher = fisher.info,fisher.invers = ginv(fisher.info),k.stu = cbind(k.aaa, k.aab, k.abb, k.bbb),k.s.tu = cbind(k.a.aa, k.a.ab, k.a.bb, k.b.aa, k.b.ab, k.b.bb)))
}

compute.augmentation <- function(components,va,vb){
  pa = ncol(va)
  pb = ncol(vb)
  n = dim(vb)[1]
  fisher = components$fisher
  k.rs = components$fisher.invers
  k.stu = components$k.stu
  k.s.tu = components$k.s.tu

  return(compute_augmentation_cpp(va, vb, fisher, k.rs, k.stu, k.s.tu))

}

compute.score <- function(x, alpha.ml, beta.ml, va, vb){

  p0p1 = getProbRD(va %*% alpha.ml, vb %*% beta.ml)
  p0 = p0p1[, 1]
  p1 = p0p1[, 2]
  n = nrow(vb)
  pA = p0
  pA[x == 1] = p1[x == 1]
  s0 = p0 * (1 - p0) # n by 1
  s1 = p1 * (1 - p1)
  sA = pA * (1 - pA)

  rho = as.vector(tanh(va * alpha.ml))

  score.alpha <- colSums(((y-pA)/sA)*(x-s0/(s0 + s1))*(1 - rho^2) * va)
  score.beta <- colSums(((y-pA)/sA)*(s0 * s1/(s0 + s1))*vb)
  return(c(score.alpha,score.beta))
}


optim.alpha <- function(alpha,beta){
  score.intial = compute.score(x,alpha,beta,va,vb)
  components = compute.components(x,alpha,beta,va,vb,weight)
  augment.intial = compute.augmentation(components,va,vb)
  return(max(abs(score.intial[1:pa] + t(augment.intial)[1:pa])))
}
optim.beta <- function(alpha,beta){
  score.intial = compute.score(x,alpha,beta,va,vb)
  components = compute.components(x,alpha,beta,va,vb,weight)
  augment.intial = compute.augmentation(components,va,vb)
  return(max(abs(score.intial[(pa+1):(pa+pb)] + t(augment.intial)[(pa+1):(pa+pb)])))
}


  Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
  alpha = alpha.start
  beta = beta.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max.step){
    step = step+1
    target.alpha <- function(alpha) { optim.alpha(alpha,beta) }
    result.a <- optim(alpha, target.alpha,control=list(maxit=max(100,max.step/10)))
    diff1 = Diff(result.a$par,alpha)
    alpha = result.a$par
    target.beta <- function(beta) { optim.beta(alpha,beta) }
    result.b <- optim(beta, target.beta,control=list(maxit=max(100,max.step/10)))
    diff  = max(diff1,Diff(result.b$par,beta))
    beta = result.b$par
  }
  opt = list(par = c(alpha,beta), convergence = (step < max.step), step = step)

  return(opt)
}



