#include <limits>
#include <vector>

#include <gluon/gluon.h>
#include <llcmcpp/go.h>

#include <fft_to_main.h>
#include <main_to_fft.h>

#include "main.decl.h"

struct main : public CBase_main, llcmcpp::Go {
  
  //TODO: cfg
  std::vector<CProxy_fft_to_main> m_from_fft;
  char** argv;
  int argc;
  uint32_t numChares;
  
  main(): FFTReady_count(0), FFTDone_count(0), printResidual_count(0) {
  }
  
  void from_fft (CProxy_fft_to_main from_fft) {
    m_from_fft.push_back(from_fft);
  }
  
  double start;
  uint64_t N;
  
  void go () {
    if (argc != 2)
      CkAbort("1 argument required\n");
    
    N = atol(argv[1]);

    if (N % numChares != 0)
      CkAbort("numChares not a factor of N\n");
    CkPrintf("Calling init\n");
    for ( int ii=0; ii<m_from_fft.size(); ++ii ) m_from_fft[ii].init(N);
  }

  uint32_t FFTReady_count;
  void FFTReady() {
    CkPrintf("FFTReady\n");
    if ( ++FFTReady_count < numChares ) return;
    FFTReady_count = 0;
    CkPrintf("FFTReady go\n");
    start = CkWallTimer();
    //TODO: Broadcast the 'go' signal to the fft chare array
    for ( int ii=0; ii<m_from_fft.size(); ++ii ) m_from_fft[ii].doFFT();
  }

  uint32_t FFTDone_count;
  void FFTDone() {
    CkPrintf("FFTDone\n");
    if ( ++FFTDone_count < numChares ) return;
    FFTDone_count = 0;
    CkPrintf("FFTDone go\n");
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("chares: %d\ncores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
             numChares, CkNumPes(), N*N, time, gflops);

    for ( int ii=0; ii<m_from_fft.size(); ++ii ) m_from_fft[ii].initValidation();
  }

  uint32_t printResidual_count;
  double printResidual_max;
  void printResidual(double r) {
    CkPrintf("printResidual\n");
    if ( printResidual_count == 0 || r > printResidual_count ) printResidual_max = r;
    if ( ++printResidual_count < numChares ) return;
    printResidual_count = 0;
    CkPrintf("residual = %g\n", printResidual_max);
    CkExit();
  }
};

#include "main.def.h"

GCMP(main)
  G_PROPERTY(int, argc);
  G_PROPERTY(char**, argv);
  G_PROPERTY(uint32_t, numChares);
  G_CHARM_PROVIDE(main_to_fft, to_fft);
  G_CPP_PROVIDE(llcmcpp::Go, go);
  G_CHARM_USE2(fft_to_main, from_fft);
GEND
