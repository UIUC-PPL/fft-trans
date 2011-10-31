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

  void FFTReady() {
    CkPrintf("FFTReady\n");
    start = CkWallTimer();
    //TODO: Broadcast the 'go' signal to the fft chare array
    for ( int ii=0; ii<m_from_fft.size(); ++ii ) m_from_fft[ii].doFFT();
  }

  void FFTDone() {
    CkPrintf("FFTDone\n");
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("chares: %d\ncores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
             numChares, CkNumPes(), N*N, time, gflops);

    for ( int ii=0; ii<m_from_fft.size(); ++ii ) m_from_fft[ii].initValidation();
  }

  void printResidual(double r) {
    CkPrintf("residual = %g\n", r);
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
