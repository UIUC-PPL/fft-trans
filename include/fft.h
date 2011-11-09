#ifndef FFT1D_FFT_H
#define FFT1D_FFT_H

#include <fftw3.h>

#include "Fft.decl.h"

struct FftCallback
{
	virtual void FFTDone ()
	{ CkPrintf("Workaround, should not get there"); CkExit(); }
};

struct Fft:
		public CBase_Fft
{
	Fft ()
	{}
	Fft ( CkMigrateMessage* )
	{}
	virtual void init ( uint64_t N, fftw_complex* input, fftw_complex* output, int direction, FftCallback* callback )
	{ CkPrintf("Workaround, should not get there"); CkExit(); }
	virtual void doFFT ()
	{ CkPrintf("Workaround, should not get there"); CkExit(); }
};

#endif // FFT1D_FFT_H
