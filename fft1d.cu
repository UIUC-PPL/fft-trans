#include "fft1d.h"

#define BLOCK_DIM 16

__global__ void twiddleKernel(complex_t* out, int N, int numChares, int k,
    double sign) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  double a, c, s, re, im;
  int idx = j*N+i;

  if (i < N && j < N/numChares) {
    a = sign * (TWOPI*(j+k*N/numChares)*i)/(N*N);
    c = cos(a);
    s = sin(a);
    re = c*out[idx].x - s*out[idx].y;
    im = s*out[idx].x + c*out[idx].y;
    out[idx].x = re;
    out[idx].y = im;
  }
}

void invokeTwiddle(complex_t* out, int N, int numChares, int k, double sign,
    cudaStream_t stream) {
  dim3 block_dim(BLOCK_DIM, BLOCK_DIM);
  dim3 grid_dim((N+block_dim.x-1) / block_dim.x,
      (N/numChares+block_dim.y-1) / block_dim.y);

  twiddleKernel<<<grid_dim, block_dim, 0, stream>>>(out, N, numChares, k, sign);
  hapiCheck(cudaPeekAtLastError());
}
