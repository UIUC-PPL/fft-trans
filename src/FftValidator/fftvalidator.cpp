#include <cstdlib>
#include <limits>

#include <fftw3.h>
#include <gluon/gluon.h>

#include <fft.h>
#include <validation.h>

#include "FftValidator.decl.h"

struct FftValidator:
		public CBase_FftValidator,
		public FftCallback
{
	//Externally configured (see bellow)
	uint32_t numChares;
	void size ( uint32_t size ) 
	{
		numChares = size;
	}
	
	//Externally configured (see bellow)
	CProxy_ValidationCallback m_callback;
	void callback ( CProxy_ValidationCallback callback )
	{
		m_callback = callback;
	}
	
	//Externally configured (see bellow)
	Fft* m_fft;
	void fft ( CProxy_Fft fft )
	{
		m_fft = fft[thisIndex].ckLocal();
	}

	FftValidator()
	{
	}
	
	FftValidator(CkMigrateMessage*)
	{
	}
	
	fftw_complex* input_buffer;
	fftw_complex* output_buffer;
	bool validating;
	uint64_t N, n;
	
	void init ( size_t N )
	{
		this->N = N;
		this->n = N*N/numChares;
		this->input_buffer = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
		this->output_buffer = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
		
		if (N % numChares != 0)
		CkAbort("numChares not a factor of N\n");
		
		srand48(thisIndex);
		for (int i = 0; i < n; i++) {
			input_buffer[i][0] = drand48();
			input_buffer[i][1] = drand48();
		}
		
		m_fft->init(N, input_buffer, output_buffer, FFTW_FORWARD, this);
		
		// Reduction to the mainchare to signal that initialization is complete
		contribute(CkCallback(CkReductionTarget(ValidationCallback, ready), m_callback));
	}
	
	void doFFT ()
	{
		validating = false;
		m_fft->doFFT();
	}
	
	void doValidation() {
		m_fft->init(N, output_buffer, output_buffer, FFTW_BACKWARD, this);
		
		validating = true;
		m_fft->doFFT();
	}
	
	void FFTDone ()
	{
		if ( !validating ) {
			// Reduction to the mainchare to signal that computation is complete
			contribute(CkCallback(CkReductionTarget(ValidationCallback, FFTDone), m_callback));
		} else {
// 			char filename[80];
// 			sprintf(filename, "%d-%ld.dump%d", numChares, N, thisIndex);
// 			writeCommFile(n, in, filename);
			
			double infNorm = 0.0;
			for (int i = 0; i < n; i++) {
				output_buffer[i][0] = output_buffer[i][0]/(N*N) - input_buffer[i][0];
				output_buffer[i][1] = output_buffer[i][1]/(N*N) - input_buffer[i][1];
				
				double mag = sqrt(pow(output_buffer[i][0], 2) + pow(output_buffer[i][1], 2));
				if (mag > infNorm) infNorm = mag;
			}
			
			double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));
			
			// Reduction to the mainchare to get the global residual
			contribute(
				sizeof(double), &r,
				CkReduction::max_double, CkCallback(CkReductionTarget(ValidationCallback, validationDone), m_callback)
			);
		}
	}
};

#include "FftValidator.def.h"

GCMP_A(FftValidator);
	G_PROPERTY2(uint32_t, size);
	G_CHARM_APROVIDE(Validation, validator);
	G_CHARM_USE2(ValidationCallback, callback);
	G_CHARM_AUSE2(Fft, fft)
GEND
