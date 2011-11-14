#include <gluon/gluon.h>
#include <llcmcpp/go.h>

#include <validation.h>

#include "ValidatorMain.decl.h"

struct ValidatorMain:
		public CBase_ValidatorMain,
		public llcmcpp::Go
{
	//Externally configured
	CProxy_Validation validator;
	
	//Externally configured
	char** argv;
	
	//Externally configured
	int argc;
	
	//Externally configured
	uint32_t numChares;
	
	ValidatorMain ()
	{
	}
	
	double start;
	uint64_t N;
	
	void go () {
		if (argc != 2)
		CkAbort("1 argument required\n");
		
		N = atol(argv[1]);

		validator.init(N);
	}

	void ready ()
	{
		start = CkWallTimer();
		validator.doFFT();
	}

	void FFTDone ()
	{
		double time = CkWallTimer() - start;
		double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
		CkPrintf("chares: %d\ncores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
				numChares, CkNumPes(), N*N, time, gflops);

		validator.doValidation();
	}

	void validationDone ( double r )
	{
		CkPrintf("residual = %g\n", r);
		CkExit();
	}
};

#include "ValidatorMain.def.h"

GCMP(ValidatorMain)
	G_PROPERTY(int, argc);
	G_PROPERTY(char**, argv);
	G_PROPERTY(uint32_t, numChares);
	G_CHARM_PROVIDE(ValidationCallback, validator_callback);
	G_CPP_PROVIDE(llcmcpp::Go, go);
	G_CHARM_AUSE(Validation, validator);
GEND
