#ifndef FFT1D_TRANSPOSE_TO_TRANSPOSE_H
#define FFT1D_TRANSPOSE_TO_TRANSPOSE_H

#include <fftw3.h>

#include "transpose_to_transpose.decl.h"

struct fftMsg : public CMessage_fftMsg {
  int source;
  fftw_complex *data;
};

struct transpose_to_transpose : public CBase_transpose_to_transpose
{
  virtual void getTranspose(fftMsg *m) = 0;
};

#endif // FFT1D_TRANSPOSE_TO_TRANSPOSE_H
