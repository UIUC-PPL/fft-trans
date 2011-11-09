#include <cstdlib>
#include <limits>

#include <fftw3.h>
#include <gluon/gluon.h>

#include <fft.h>
#include <transpose.h>

#include "FftComputer.decl.h"

#define TWOPI 6.283185307179586

struct FftComputer:
		public CBase_FftComputer,
		public transpose_callback
{
	FftComputer_SDAG_CODE
	
	//Externally configured (see bellow)
	uint32_t numChares;
	void size ( uint32_t size ) 
	{
		numChares = size;
	}
	
	//Externally configured (see bellow)
	transpose* m_transpose;
	void transposer ( CProxy_transpose transposer )
	{
		m_transpose = transposer[thisIndex].ckLocal();
	}
	
	
	FftComputer ():
			p1(0)
	{
		__sdag_init();
	}
	
	FftComputer ( CkMigrateMessage* ):
			p1(0)
	{
		__sdag_init();
	}
	
	FftCallback* callback;
	int iteration, count;
	fftw_plan p1;
	fftw_complex *in, *out;
	uint64_t N;
	double sign;
	
	void init ( uint64_t N, fftw_complex* input, fftw_complex* output, int direction, FftCallback* callback )
	{
		m_transpose->init(N);
		
		if (p1) fftw_destroy_plan(p1);
		
		this->N = N;
		this->in = input;
		this->out = output;
		this->sign = direction;
		this->callback = callback;
		
		int length[] = {N};
		this->p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
		                              out, length, 1, N, direction, FFTW_ESTIMATE);
	}
	
	void twiddle()
	{
		double a, c, s, re, im;
		
		int k = thisIndex;
		for (int i = 0; i < N/numChares; i++)
		for (int j = 0; j < N; j++) {
			a = sign * (TWOPI*(i+k*N/numChares)*j)/(N*N);
			c = cos(a);
			s = sin(a);
			
			int idx = i*N+j;
			
			re = c*out[idx][0] - s*out[idx][1];
			im = s*out[idx][0] + c*out[idx][1];
			out[idx][0] = re;
			out[idx][1] = im;
		}
	}
};

#include "FftComputer.def.h"

GCMP_A(FftComputer);
	G_PROPERTY2(uint32_t, size);
	G_CHARM_APROVIDE(Fft, fft);
	G_CHARM_AUSE2(transpose, transposer)
GEND
