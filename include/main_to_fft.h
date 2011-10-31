#ifndef FFT1D_MAIN_TO_FFT_H
#define FFT1D_MAIN_TO_FFT_H

#include "main_to_fft.decl.h"

struct main_to_fft : public CBase_main_to_fft
{
  virtual void FFTReady() = 0;
  virtual void FFTDone() = 0;
  virtual void printResidual(double residual) = 0;
};

#endif // FFT1D_MAIN_TO_FFT_H
