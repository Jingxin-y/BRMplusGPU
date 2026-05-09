// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @importFrom Rcpp evalCpp
//' @useDynLib brm
//' @exportPattern ˆ[[:alpha:]]+
// [[Rcpp::export]]
arma::vec compute_augmentation_cpp(
  const arma::mat& va,
  const arma::mat& vb,
  const arma::mat& fisher,
  const arma::mat& k_rs,
  const arma::mat& k_stu,
  const arma::mat& k_s_tu) {

  int pa = static_cast<int>(va.n_cols);
  int pb = static_cast<int>(vb.n_cols);
  int n  = static_cast<int>(vb.n_rows);

  arma::mat kaa(n, pa, arma::fill::zeros);
  arma::mat kab(n, pa, arma::fill::zeros);
  arma::mat kba(n, pb, arma::fill::zeros);
  arma::mat kbb(n, pb, arma::fill::zeros);
  arma::mat b1_a(n, pa, arma::fill::zeros);
  arma::mat b1_b(n, pb, arma::fill::zeros);

  // a-block
  for (int a1 = 0; a1 < pa; ++a1) {
    kaa.col(a1).zeros();
    for (int a2 = 0; a2 < pa; ++a2) {
      arma::vec kaa_m(n, arma::fill::zeros);

      for (int a3 = 0; a3 < pa; ++a3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kaa_m += k_rs(a3, a4) *
            (k_stu.col(0) + k_s_tu.col(0)) %
            va.col(a2) % va.col(a3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kaa_m += k_rs(a3, 1 + b4) *
            (k_stu.col(1) + k_s_tu.col(1)) %
            va.col(a2) % va.col(a3) % vb.col(b4);
        }
      }
      for (int b3 = 0; b3 < pb; ++b3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kaa_m += k_rs(1 + b3, a4) *
            (k_stu.col(1) + k_s_tu.col(1)) %
            va.col(a2) % vb.col(b3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kaa_m += k_rs(1 + b3, 1 + b4) *
            (k_stu.col(2) + k_s_tu.col(2)) %
            va.col(a2) % vb.col(b3) % vb.col(b4);
        }
      }
      kaa.col(a1) += k_rs(a1, a2) * kaa_m;
    }

    kab.col(a1).zeros();
    for (int b2 = 0; b2 < pb; ++b2) {
      arma::vec kab_m(n, arma::fill::zeros);

      for (int a3 = 0; a3 < pa; ++a3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kab_m += k_rs(a3, a4) *
            (k_stu.col(1) + k_s_tu.col(3)) %
            vb.col(b2) % va.col(a3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kab_m += k_rs(a3, 1 + b4) *
            (k_stu.col(2) + k_s_tu.col(4)) %
            vb.col(b2) % va.col(a3) % vb.col(b4);
        }
      }
      for (int b3 = 0; b3 < pb; ++b3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kab_m += k_rs(1 + b3, a4) *
            (k_stu.col(2) + k_s_tu.col(4)) %
            vb.col(b2) % vb.col(b3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kab_m += k_rs(1 + b3, 1 + b4) *
            (k_stu.col(3) + k_s_tu.col(5)) %
            vb.col(b2) % vb.col(b3) % vb.col(b4);
        }
      }
      kab.col(a1) += k_rs(a1, 1 + b2) * kab_m;
    }

    b1_a.col(a1) = kaa.col(a1) + kab.col(a1);
  }

  // b-block
  for (int b1 = 0; b1 < pb; ++b1) {
    kba.col(b1).zeros();
    for (int a2 = 0; a2 < pa; ++a2) {
      arma::vec kba_m(n, arma::fill::zeros);

      for (int a3 = 0; a3 < pa; ++a3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kba_m += k_rs(a3, a4) *
            (k_stu.col(0) + k_s_tu.col(0)) %
            va.col(a2) % va.col(a3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kba_m += k_rs(a3, 1 + b4) *
            (k_stu.col(1) + k_s_tu.col(1)) %
            va.col(a2) % va.col(a3) % vb.col(b4);
        }
      }
      for (int b3 = 0; b3 < pb; ++b3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kba_m += k_rs(1 + b3, a4) *
            (k_stu.col(1) + k_s_tu.col(1)) %
            va.col(a2) % vb.col(b3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kba_m += k_rs(1 + b3, b4) *
            (k_stu.col(2) + k_s_tu.col(2)) %
            va.col(a2) % vb.col(b3) % vb.col(b4);
        }
      }
      kba.col(b1) += k_rs(1 + b1, a2) * kba_m;
    }

    kbb.col(b1).zeros();
    for (int b2 = 0; b2 < pb; ++b2) {
      arma::vec kbb_m(n, arma::fill::zeros);

      for (int a3 = 0; a3 < pa; ++a3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kbb_m += k_rs(a3, a4) *
            (k_stu.col(1) + k_s_tu.col(3)) %
            vb.col(b2) % va.col(a3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kbb_m += k_rs(a3, 1 + b4) *
            (k_stu.col(2) + k_s_tu.col(4)) %
            vb.col(b2) % va.col(a3) % vb.col(b4);
        }
      }
      for (int b3 = 0; b3 < pb; ++b3) {
        for (int a4 = 0; a4 < pa; ++a4) {
          kbb_m += k_rs(1 + b3, a4) *
            (k_stu.col(2) + k_s_tu.col(4)) %
            vb.col(b2) % vb.col(b3) % va.col(a4);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kbb_m = kbb_m + k_rs(1 + b3, 1 + b4) *
            (k_stu.col(3) + k_s_tu.col(5)) %
            vb.col(b2) % vb.col(b3) % vb.col(b4);
        }
      }
      kbb.col(b1) += k_rs(1 + b1, 1 + b2) * kbb_m;
    }

    b1_b.col(b1) = kba.col(b1) + kbb.col(b1);
  }

  // colMeans
  arma::rowvec b1_a_mean = arma::mean(b1_a, 0);
  arma::rowvec b1_b_mean = arma::mean(b1_b, 0);

  // b1 = -c(b1.a, b1.b)/2
  arma::vec b1 = -0.5 * join_cols(b1_a_mean.t(), b1_b_mean.t());

  // expect.A = -fisher %*% b1 / n
  arma::vec expect_A = -(fisher * b1) / static_cast<double>(n);

  return expect_A;
}
