#ifndef __FFT1D_H_
#define __FFT1D_H_

// Set CPU as default execution mode
#if !defined MODE_CUDA
#define MODE_CPU
#endif

#ifdef MODE_CPU
#include <fftw3.h>
typedef fftw_complex complex_t;
#elif defined MODE_CUDA
#include <cufft.h>
#include <cuda_runtime_api.h>
#include "hapi.h"
typedef cufftComplex complex_t;
#endif

#endif // __FFT1D_H_
