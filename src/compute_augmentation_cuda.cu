#include <cuda_runtime.h>

#include <cmath>
#include <cstdlib>
#include <vector>

namespace {

inline int k_index(int row, int col, int nrow) {
  return row + col * nrow;
}

__device__ double mat_get(const double* x, int row, int col, int nrow) {
  return x[row + col * nrow];
}

__global__ void augmentation_rows_kernel(
  const double* va,
  const double* vb,
  const double* k_rs,
  const double* k_stu,
  const double* k_s_tu,
  int n,
  int pa,
  int pb,
  double* rows) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int out_col = blockIdx.y;
  int kdim = pa + pb;

  if (i >= n) {
    return;
  }

  double value = 0.0;

  if (out_col < pa) {
    int a1 = out_col;
    double kaa = 0.0;
    double kab = 0.0;

    for (int a2 = 0; a2 < pa; ++a2) {
      double kaa_m = 0.0;
      double va_a2 = mat_get(va, i, a2, n);

      for (int a3 = 0; a3 < pa; ++a3) {
        double va_a3 = mat_get(va, i, a3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kaa_m += mat_get(k_rs, a3, a4, kdim) *
            (mat_get(k_stu, i, 0, n) + mat_get(k_s_tu, i, 0, n)) *
            va_a2 * va_a3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kaa_m += mat_get(k_rs, a3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 1, n) + mat_get(k_s_tu, i, 1, n)) *
            va_a2 * va_a3 * mat_get(vb, i, b4, n);
        }
      }

      for (int b3 = 0; b3 < pb; ++b3) {
        double vb_b3 = mat_get(vb, i, b3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kaa_m += mat_get(k_rs, 1 + b3, a4, kdim) *
            (mat_get(k_stu, i, 1, n) + mat_get(k_s_tu, i, 1, n)) *
            va_a2 * vb_b3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kaa_m += mat_get(k_rs, 1 + b3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 2, n) + mat_get(k_s_tu, i, 2, n)) *
            va_a2 * vb_b3 * mat_get(vb, i, b4, n);
        }
      }
      kaa += mat_get(k_rs, a1, a2, kdim) * kaa_m;
    }

    for (int b2 = 0; b2 < pb; ++b2) {
      double kab_m = 0.0;
      double vb_b2 = mat_get(vb, i, b2, n);

      for (int a3 = 0; a3 < pa; ++a3) {
        double va_a3 = mat_get(va, i, a3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kab_m += mat_get(k_rs, a3, a4, kdim) *
            (mat_get(k_stu, i, 1, n) + mat_get(k_s_tu, i, 3, n)) *
            vb_b2 * va_a3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kab_m += mat_get(k_rs, a3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 2, n) + mat_get(k_s_tu, i, 4, n)) *
            vb_b2 * va_a3 * mat_get(vb, i, b4, n);
        }
      }

      for (int b3 = 0; b3 < pb; ++b3) {
        double vb_b3 = mat_get(vb, i, b3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kab_m += mat_get(k_rs, 1 + b3, a4, kdim) *
            (mat_get(k_stu, i, 2, n) + mat_get(k_s_tu, i, 4, n)) *
            vb_b2 * vb_b3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kab_m += mat_get(k_rs, 1 + b3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 3, n) + mat_get(k_s_tu, i, 5, n)) *
            vb_b2 * vb_b3 * mat_get(vb, i, b4, n);
        }
      }
      kab += mat_get(k_rs, a1, 1 + b2, kdim) * kab_m;
    }

    value = kaa + kab;
  } else {
    int b1 = out_col - pa;
    double kba = 0.0;
    double kbb = 0.0;

    for (int a2 = 0; a2 < pa; ++a2) {
      double kba_m = 0.0;
      double va_a2 = mat_get(va, i, a2, n);

      for (int a3 = 0; a3 < pa; ++a3) {
        double va_a3 = mat_get(va, i, a3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kba_m += mat_get(k_rs, a3, a4, kdim) *
            (mat_get(k_stu, i, 0, n) + mat_get(k_s_tu, i, 0, n)) *
            va_a2 * va_a3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kba_m += mat_get(k_rs, a3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 1, n) + mat_get(k_s_tu, i, 1, n)) *
            va_a2 * va_a3 * mat_get(vb, i, b4, n);
        }
      }

      for (int b3 = 0; b3 < pb; ++b3) {
        double vb_b3 = mat_get(vb, i, b3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kba_m += mat_get(k_rs, 1 + b3, a4, kdim) *
            (mat_get(k_stu, i, 1, n) + mat_get(k_s_tu, i, 1, n)) *
            va_a2 * vb_b3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kba_m += mat_get(k_rs, 1 + b3, b4, kdim) *
            (mat_get(k_stu, i, 2, n) + mat_get(k_s_tu, i, 2, n)) *
            va_a2 * vb_b3 * mat_get(vb, i, b4, n);
        }
      }
      kba += mat_get(k_rs, 1 + b1, a2, kdim) * kba_m;
    }

    for (int b2 = 0; b2 < pb; ++b2) {
      double kbb_m = 0.0;
      double vb_b2 = mat_get(vb, i, b2, n);

      for (int a3 = 0; a3 < pa; ++a3) {
        double va_a3 = mat_get(va, i, a3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kbb_m += mat_get(k_rs, a3, a4, kdim) *
            (mat_get(k_stu, i, 1, n) + mat_get(k_s_tu, i, 3, n)) *
            vb_b2 * va_a3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kbb_m += mat_get(k_rs, a3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 2, n) + mat_get(k_s_tu, i, 4, n)) *
            vb_b2 * va_a3 * mat_get(vb, i, b4, n);
        }
      }

      for (int b3 = 0; b3 < pb; ++b3) {
        double vb_b3 = mat_get(vb, i, b3, n);
        for (int a4 = 0; a4 < pa; ++a4) {
          kbb_m += mat_get(k_rs, 1 + b3, a4, kdim) *
            (mat_get(k_stu, i, 2, n) + mat_get(k_s_tu, i, 4, n)) *
            vb_b2 * vb_b3 * mat_get(va, i, a4, n);
        }
        for (int b4 = 0; b4 < pb; ++b4) {
          kbb_m += mat_get(k_rs, 1 + b3, 1 + b4, kdim) *
            (mat_get(k_stu, i, 3, n) + mat_get(k_s_tu, i, 5, n)) *
            vb_b2 * vb_b3 * mat_get(vb, i, b4, n);
        }
      }
      kbb += mat_get(k_rs, 1 + b1, 1 + b2, kdim) * kbb_m;
    }

    value = kba + kbb;
  }

  rows[i + out_col * n] = value;
}

int fail(cudaError_t status) {
  return status == cudaSuccess ? 0 : static_cast<int>(status);
}

} // namespace

extern "C" int compute_augmentation_cuda(
  const double* va,
  const double* vb,
  const double* fisher,
  const double* k_rs,
  const double* k_stu,
  const double* k_s_tu,
  int n,
  int pa,
  int pb,
  double* out) {

  int kdim = pa + pb;
  std::size_t va_bytes = static_cast<std::size_t>(n) * pa * sizeof(double);
  std::size_t vb_bytes = static_cast<std::size_t>(n) * pb * sizeof(double);
  std::size_t kdim_bytes = static_cast<std::size_t>(kdim) * kdim * sizeof(double);
  std::size_t fisher_bytes = kdim_bytes;
  std::size_t k_stu_bytes = static_cast<std::size_t>(n) * 4 * sizeof(double);
  std::size_t k_s_tu_bytes = static_cast<std::size_t>(n) * 6 * sizeof(double);
  std::size_t rows_bytes = static_cast<std::size_t>(n) * kdim * sizeof(double);

  double* d_va = nullptr;
  double* d_vb = nullptr;
  double* d_k_rs = nullptr;
  double* d_k_stu = nullptr;
  double* d_k_s_tu = nullptr;
  double* d_rows = nullptr;

  if (cudaMalloc(&d_va, va_bytes) != cudaSuccess ||
      cudaMalloc(&d_vb, vb_bytes) != cudaSuccess ||
      cudaMalloc(&d_k_rs, kdim_bytes) != cudaSuccess ||
      cudaMalloc(&d_k_stu, k_stu_bytes) != cudaSuccess ||
      cudaMalloc(&d_k_s_tu, k_s_tu_bytes) != cudaSuccess ||
      cudaMalloc(&d_rows, rows_bytes) != cudaSuccess) {
    cudaFree(d_va);
    cudaFree(d_vb);
    cudaFree(d_k_rs);
    cudaFree(d_k_stu);
    cudaFree(d_k_s_tu);
    cudaFree(d_rows);
    return 1;
  }

  int status = 0;
  status |= fail(cudaMemcpy(d_va, va, va_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_vb, vb, vb_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_k_rs, k_rs, kdim_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_k_stu, k_stu, k_stu_bytes, cudaMemcpyHostToDevice));
  status |= fail(cudaMemcpy(d_k_s_tu, k_s_tu, k_s_tu_bytes, cudaMemcpyHostToDevice));

  if (status == 0) {
    dim3 block(256);
    dim3 grid((n + block.x - 1) / block.x, kdim);
    augmentation_rows_kernel<<<grid, block>>>(d_va, d_vb, d_k_rs, d_k_stu, d_k_s_tu, n, pa, pb, d_rows);
    status |= fail(cudaGetLastError());
    status |= fail(cudaDeviceSynchronize());
  }

  std::vector<double> rows(static_cast<std::size_t>(n) * kdim);
  if (status == 0) {
    status |= fail(cudaMemcpy(rows.data(), d_rows, rows_bytes, cudaMemcpyDeviceToHost));
  }

  cudaFree(d_va);
  cudaFree(d_vb);
  cudaFree(d_k_rs);
  cudaFree(d_k_stu);
  cudaFree(d_k_s_tu);
  cudaFree(d_rows);

  if (status != 0) {
    return status;
  }

  std::vector<double> b1(kdim, 0.0);
  for (int col = 0; col < kdim; ++col) {
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      total += rows[k_index(i, col, n)];
    }
    b1[col] = -0.5 * total / static_cast<double>(n);
  }

  for (int row = 0; row < kdim; ++row) {
    double total = 0.0;
    for (int col = 0; col < kdim; ++col) {
      total += fisher[k_index(row, col, kdim)] * b1[col];
    }
    out[row] = -total / static_cast<double>(n);
  }

  return 0;
}
