#ifndef FFT1D_FFT_TO_MAIN_H
#define FFT1D_FFT_TO_MAIN_H

#include "fft_to_main.decl.h"

struct fft_to_main : public CBase_fft_to_main
{
  virtual void init ( uint64_t N ) = 0;
  virtual void doFFT () = 0;
  virtual void initValidation() = 0;
};

#endif // FFT1D_FFT_TO_MAIN_H
