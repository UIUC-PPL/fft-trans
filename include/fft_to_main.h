#ifndef FFT1D_FFT_TO_MAIN_H
#define FFT1D_FFT_TO_MAIN_H

#include "fft_to_main.decl.h"

struct fft_to_main : public CBase_fft_to_main
{
  fft_to_main() {}
  fft_to_main(CkMigrateMessage*) { CkPrintf("Ugly Workaround, should not get there"); CkExit(); }
  virtual void init ( uint64_t N, uint32_t numChares ) { CkPrintf("Ugly Workaround, should not get there"); CkExit(); }
  virtual void doFFT () { CkPrintf("Ugly Workaround, should not get there"); CkExit(); }
  virtual void initValidation() { CkPrintf("Ugly Workaround, should not get there"); CkExit(); }
};

#endif // FFT1D_FFT_TO_MAIN_H
