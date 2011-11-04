#ifndef FFT1D_TRANSPOSE_TO_FFT_H
#define FFT1D_TRANSPOSE_TO_FFT_H

#include <fftw3.h>

struct fft_to_transpose
{
  virtual void transposeDone (int) = 0;
};

struct transpose_to_fft
{
  virtual void init(uint64_t N) = 0;
  virtual void sendTranspose(int iteration, fftw_complex* in, fftw_complex* out, fft_to_transpose& callback) = 0;
};

#endif // FFT1D_TRANSPOSE_TO_FFT_H
