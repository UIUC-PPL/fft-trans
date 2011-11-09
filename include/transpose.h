#ifndef FFT1D_TRANSPOSE_TO_FFT_H
#define FFT1D_TRANSPOSE_TO_FFT_H

#include <stdint.h>

#include <fftw3.h>

#include "transpose.decl.h"

struct transpose_callback
{
	virtual void transposeDone (int) = 0;
};

struct transpose: public CBase_transpose
{
	virtual void init(uint64_t N) = 0;
	virtual void sendTranspose(int iteration, fftw_complex* input, fftw_complex* output, transpose_callback* cb) = 0;
};

#endif // FFT1D_TRANSPOSE_TO_FFT_H
