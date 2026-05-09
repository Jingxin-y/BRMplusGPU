#include "RcppArmadillo.h"
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

#ifdef BRM_USE_CUDA
extern "C" int exact_nll_cuda(
  int is_rr,
  const double* y,
  const double* x,
  const double* va,
  const double* vb,
  const double* alpha,
  const double* beta,
  const double* weight,
  int n,
  int pa,
  int pb,
  double eps,
  double* out);
#endif

namespace {

bool use_gpu_requested_for_exact() {
  const char* value = std::getenv("BRM_USE_GPU");
  if (value == nullptr) {
    return false;
  }

  std::string flag(value);
  return flag == "1" || flag == "true" || flag == "TRUE" || flag == "yes" || flag == "YES";
}

double get_prb_aux(double x) {
  if (x < 17.0 && x > -500.0) {
    return 0.5 * std::exp(x) * (-1.0 + std::sqrt(1.0 + 4.0 * std::exp(-x)));
  }
  return x < 0.0 ? 0.0 : 1.0;
}

bool same(double x, double y) {
  return std::abs(x - y) < std::sqrt(std::numeric_limits<double>::epsilon());
}

void get_prob_rr(double logrr, double logop, double& p0, double& p1) {
  if (logop < -12.0 || logop > 12.0 || logrr < -12.0 || logrr > 12.0) {
    if (logrr < -12.0 || (logop < -12.0 && logrr < 0.0)) {
      p0 = get_prb_aux(logop - logrr);
      p1 = 0.0;
    } else if (logrr > 12.0 || (logop < -12.0 && logrr > 0.0)) {
      p0 = 0.0;
      p1 = get_prb_aux(logop + logrr);
    } else {
      p0 = std::min(std::exp(-logrr), 1.0);
      p1 = std::min(std::exp(logrr), 1.0);
    }
  } else {
    if (same(logop, 0.0)) {
      p0 = 1.0 / (1.0 + std::exp(logrr));
    } else {
      double exp_logrr = std::exp(logrr);
      double exp_logop = std::exp(logop);
      p0 = (-(exp_logrr + 1.0) * exp_logop +
        std::sqrt(std::exp(2.0 * logop) * std::pow(exp_logrr + 1.0, 2.0) +
          4.0 * std::exp(logrr + logop) * (1.0 - exp_logop))) /
        (2.0 * exp_logrr * (1.0 - exp_logop));
    }
    p1 = std::exp(logrr) * p0;
  }
}

void get_prob_rd(double atanhrd, double logop, double& p0, double& p1) {
  double rd = std::tanh(atanhrd);
  if (logop > 12.0) {
    if (atanhrd < 0.0) {
      p0 = 1.0;
      p1 = p0 + rd;
    } else {
      p1 = 1.0;
      p0 = p1 - rd;
    }
  } else if (logop < -12.0) {
    if (atanhrd < 0.0) {
      p0 = -rd;
      p1 = 0.0;
    } else {
      p1 = rd;
      p0 = 0.0;
    }
  } else {
    if (same(logop, 0.0)) {
      p0 = 0.5 * (1.0 - rd);
    } else {
      double exp_logop = std::exp(logop);
      double tmp = exp_logop * (rd - 2.0) - rd;
      p0 = (-tmp - std::sqrt(tmp * tmp + 4.0 * exp_logop * (1.0 - rd) * (1.0 - exp_logop))) /
        (2.0 * (exp_logop - 1.0));
    }
    p1 = p0 + rd;
  }
}

double exact_nll_cpu(
  bool is_rr,
  const arma::vec& y,
  const arma::vec& x,
  const arma::mat& va,
  const arma::mat& vb,
  const arma::vec& alpha,
  const arma::vec& beta,
  const arma::vec& weight,
  double eps) {

  double total = 0.0;
  int n = static_cast<int>(y.n_elem);
  int pa = static_cast<int>(alpha.n_elem);
  int pb = static_cast<int>(beta.n_elem);

  for (int i = 0; i < n; ++i) {
    double eta_a = 0.0;
    double eta_b = 0.0;

    for (int j = 0; j < pa; ++j) {
      eta_a += va(i, j) * alpha[j];
    }
    for (int j = 0; j < pb; ++j) {
      eta_b += vb(i, j) * beta[j];
    }

    double p0 = 0.0;
    double p1 = 0.0;
    if (is_rr) {
      get_prob_rr(eta_a, eta_b, p0, p1);
    } else {
      get_prob_rd(eta_a, eta_b, p0, p1);
    }

    double p = x[i] == 0.0 ? p0 : p1;
    p = std::min(std::max(p, eps), 1.0 - eps);
    total -= ((1.0 - y[i]) * std::log(1.0 - p) + y[i] * std::log(p)) * weight[i];
  }

  return total;
}

} // namespace

//' @importFrom Rcpp evalCpp
//' @useDynLib brm
//' @exportPattern 藛[[:alpha:]]+
// [[Rcpp::export]]
double exact_nll_cpp(
  std::string param,
  const arma::vec& y,
  const arma::vec& x,
  const arma::mat& va,
  const arma::mat& vb,
  const arma::vec& alpha,
  const arma::vec& beta,
  const arma::vec& weight,
  double eps = 1e-12) {

  bool is_rr = param == "RR";

  if (!use_gpu_requested_for_exact()) {
    return exact_nll_cpu(is_rr, y, x, va, vb, alpha, beta, weight, eps);
  }

#ifdef BRM_USE_CUDA
  double out = 0.0;
  int status = exact_nll_cuda(
    is_rr ? 1 : 0,
    y.memptr(),
    x.memptr(),
    va.memptr(),
    vb.memptr(),
    alpha.memptr(),
    beta.memptr(),
    weight.memptr(),
    static_cast<int>(y.n_elem),
    static_cast<int>(alpha.n_elem),
    static_cast<int>(beta.n_elem),
    eps,
    &out);

  if (status == 0 && std::isfinite(out)) {
    return out;
  }

  static bool warned_cuda_failure = false;
  if (!warned_cuda_failure) {
    Rcpp::warning("CUDA exact likelihood failed; falling back to CPU implementation.");
    warned_cuda_failure = true;
  }
#else
  static bool warned_no_cuda = false;
  if (!warned_no_cuda) {
    Rcpp::warning("BRM_USE_GPU is set, but brm was not compiled with CUDA; falling back to CPU exact likelihood.");
    warned_no_cuda = true;
  }
#endif

  return exact_nll_cpu(is_rr, y, x, va, vb, alpha, beta, weight, eps);
}
