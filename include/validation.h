#ifndef FFT1D_VALIDATION_H
#define FFT1D_VALIDATION_H

#include "Validation.decl.h"

struct ValidationCallback : public CBase_ValidationCallback
{
	virtual void ready () = 0;
	virtual void FFTDone () = 0;
	virtual void validationDone ( double residual ) = 0;
};

struct Validation : public CBase_Validation
{
	Validation ()
	{}
	
	Validation ( CkMigrateMessage* )
	{}
	
	virtual void init ( uint64_t N )
	{ CkPrintf("Workaround implementation, should not get there"); CkExit(); }
	
	virtual void doFFT ()
	{ CkPrintf("Workaround implementation, should not get there"); CkExit(); }
	
	virtual void doValidation ()
	{ CkPrintf("Workaround implementation, should not get there"); CkExit(); }
	
};

#endif // FFT1D_VALIDATION_H
