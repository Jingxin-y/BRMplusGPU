#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace {

__device__ double mat_get(const double* x, int row, int col, int nrow) {
  return x[row + col * nrow];
}

__device__ bool same_device(double x, double y) {
  return fabs(x - y) < 1.4901161193847656e-8;
}

__device__ double get_prb_aux_device(double x) {
  if (x < 17.0 && x > -500.0) {
    return 0.5 * exp(x) * (-1.0 + sqrt(1.0 + 4.0 * exp(-x)));
  }
  return x < 0.0 ? 0.0 : 1.0;
}

__device__ void get_prob_rr_device(double logrr, double logop, double* p0, double* p1) {
  if (logop < -12.0 || logop > 12.0 || logrr < -12.0 || logrr > 12.0) {
    if (logrr < -12.0 || (logop < -12.0 && logrr < 0.0)) {
      *p0 = get_prb_aux_device(logop - logrr);
      *p1 = 0.0;
    } else if (logrr > 12.0 || (logop < -12.0 && logrr > 0.0)) {
      *p0 = 0.0;
      *p1 = get_prb_aux_device(logop + logrr);
    } else {
      *p0 = fmin(exp(-logrr), 1.0);
      *p1 = fmin(exp(logrr), 1.0);
    }
  } else {
    if (same_device(logop, 0.0)) {
      *p0 = 1.0 / (1.0 + exp(logrr));
    } else {
      double exp_logrr = exp(logrr);
      double exp_logop = exp(logop);
      *p0 = (-(exp_logrr + 1.0) * exp_logop +
        sqrt(exp(2.0 * logop) * pow(exp_logrr + 1.0, 2.0) +
          4.0 * exp(logrr + logop) * (1.0 - exp_logop))) /
        (2.0 * exp_logrr * (1.0 - exp_logop));
    }
    *p1 = exp(logrr) * (*p0);
  }
}

__device__ void get_prob_rd_device(double atanhrd, double logop, double* p0, double* p1) {
  double rd = tanh(atanhrd);
  if (logop > 12.0) {
    if (atanhrd < 0.0) {
      *p0 = 1.0;
      *p1 = *p0 + rd;
    } else {
      *p1 = 1.0;
      *p0 = *p1 - rd;
    }
  } else if (logop < -12.0) {
    if (atanhrd < 0.0) {
      *p0 = -rd;
      *p1 = 0.0;
    } else {
      *p1 = rd;
      *p0 = 0.0;
    }
  } else {
    if (same_device(logop, 0.0)) {
      *p0 = 0.5 * (1.0 - rd);
    } else {
      double exp_logop = exp(logop);
      double tmp = exp_logop * (rd - 2.0) - rd;
      *p0 = (-tmp - sqrt(tmp * tmp + 4.0 * exp_logop * (1.0 - rd) * (1.0 - exp_logop))) /
        (2.0 * (exp_logop - 1.0));
    }
    *p1 = *p0 + rd;
  }
}

__global__ void exact_nll_kernel(
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
  double* contributions) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) {
    return;
  }

  double eta_a = 0.0;
  double eta_b = 0.0;
  for (int j = 0; j < pa; ++j) {
    eta_a += mat_get(va, i, j, n) * alpha[j];
  }
  for (int j = 0; j < pb; ++j) {
    eta_b += mat_get(vb, i, j, n) * beta[j];
  }

  double p0 = 0.0;
  double p1 = 0.0;
  if (is_rr == 1) {
    get_prob_rr_device(eta_a, eta_b, &p0, &p1);
  } else {
    get_prob_rd_device(eta_a, eta_b, &p0, &p1);
  }

  double p = x[i] == 0.0 ? p0 : p1;
  p = fmin(fmax(p, eps), 1.0 - eps);
  contributions[i] = -((1.0 - y[i]) * log(1.0 - p) + y[i] * log(p)) * weight[i];
}

int fail(cudaError_t status) {
  return status == cudaSuccess ? 0 : static_cast<int>(status);
}

} // namespace

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
  double* out) {

  std::size_t n_bytes = static_cast<std::size_t>(n) * sizeof(double);
  std::size_t va_bytes = static_cast<std::size_t>(n) * pa * sizeof(double);
  std::size_t vb_bytes = static_cast<std::size_t>(n) * pb * sizeof(double);
  std::size_t alpha_bytes = static_cast<std::size_t>(pa) * sizeof(double);
  std::size_t beta_bytes = static_cast<std::size_t>(pb) * sizeof(double);

  double* d_y = nullptr;
  double* d_x = nullptr;
  double* d_va = nullptr;
  double* d_vb = nullptr;
  double* d_alpha = nullptr;
  double* d_beta = nullptr;
  double* d_weight = nullptr;
  double* d_contrib = nullptr;

  if (cudaMalloc(&d_y, n_bytes) != cudaSuccess ||
      cudaMalloc(&d_x, n_bytes) != cudaSuccess ||
      cudaMalloc(&d_va, va_bytes) != cudaSuccess ||
      cudaMalloc(&d_vb, vb_bytes) != cudaSuccess ||
      cudaMalloc(&d_alpha, alpha_bytes) != cudaSuccess ||
      cudaMalloc(&d_beta, beta_bytes) != cudaSuccess ||
      cudaMalloc(&d_weight, n_bytes) != cudaSuccess ||
      cudaMalloc(&d_contrib, n_bytes) != cudaSuccess) {
    cudaFree(d_y);
    cudaFree(d_x);
    cudaFree(d_va);
    cudaFree(d_vb);
    cudaFree(d_alpha);
    cudaFree(d_beta);
    cudaFree(d_weight);
    cudaFree(d_contrib);
    return 1;
  }

  int status = 0;
  status |= fail(cudaMemcpy(d_y, y, n_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_x, x, n_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_va, va, va_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_vb, vb, vb_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_alpha, alpha, alpha_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_beta, beta, beta_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_weight, weight, n_bytes, cudaMemcpyHostToDevice));

  if (status == 0) {
    dim3 block(256);
    dim3 grid((n + block.x - 1) / block.x);
    exact_nll_kernel<<<grid, block>>>(is_rr, d_y, d_x, d_va, d_vb, d_alpha, d_beta, d_weight, n, pa, pb, eps, d_contrib);
    status |= fail(cudaGetLastError());
    status |= fail(cudaDeviceSynchronize());
  }

  std::vector<double> contributions(n);
  if (status == 0) {
    status |= fail(cudaMemcpy(contributions.data(), d_contrib, n_bytes, cudaMemcpyDeviceToHost));
  }

  cudaFree(d_y);
  cudaFree(d_x);
  cudaFree(d_va);
  cudaFree(d_vb);
  cudaFree(d_alpha);
  cudaFree(d_beta);
  cudaFree(d_weight);
  cudaFree(d_contrib);

  if (status != 0) {
    return status;
  }

  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    total += contributions[i];
  }

  *out = total;
  return 0;
}
