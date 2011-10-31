#ifndef FFT1D_TRANSPOSE_TO_FFT_H
#define FFT1D_TRANSPOSE_TO_FFT_H

#include <fftw3.h>

struct transpose_to_fft
{
  virtual void init(uint64_t N) = 0;
  virtual void sendTranspose(int iteration, fftw_complex* in, fftw_complex* out) = 0;
};

#endif // FFT1D_TRANSPOSE_TO_FFT_H
