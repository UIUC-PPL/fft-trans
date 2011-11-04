#include <cstdlib>
#include <limits>

#include <fftw3.h>
#include <gluon/gluon.h>

#include <fft_to_main.h>
#include <fileio.h>
#include <main_to_fft.h>
#include <transpose.h>

#include "fft.decl.h"

#define TWOPI 6.283185307179586

struct fft : public CBase_fft, transpose_callback {
  fft(CkMigrateMessage*) { CkPrintf("Ugly Workaround, should not get there"); CkExit(); }
  fft_SDAG_CODE
  
  //TODO: cfg
  CProxy_main_to_fft m_main;
  transpose* m_transpose;
  uint32_t numChares;
  
  int iteration, count;
  fftw_plan p1;
  fftw_complex *in, *out;
  uint64_t N, n;
  bool validating;
  
  fft() {
    __sdag_init();
  }
  
  void from_transpose ( CProxy_transpose from_transpose )
  {
    m_transpose = from_transpose[thisIndex].ckLocal();
  }
  
  void from_main ( CProxy_main_to_fft from_main )
  {
    m_main = from_main;
  }
  
  void init(uint64_t N) {
    m_transpose->init(N);
    
    validating = false;
    
    this->N = N;
    n = N*N/numChares;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    
    int length[] = {N};
    p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
                            out, length, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);
    
    srand48(thisIndex);
    for (int i = 0; i < n; i++) {
      in[i][0] = drand48();
      in[i][1] = drand48();
    }
    
    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(main_to_fft,FFTReady), m_main));
    
  }
  
  void twiddle(double sign) {
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
  
  void initValidation() {
    memcpy(in, out, sizeof(fftw_complex) * n);
    
    validating = true;
    fftw_destroy_plan(p1);
    int length[] = {N};
    p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
                            out, length, 1, N, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    contribute(CkCallback(CkReductionTarget(main_to_fft,FFTReady), m_main));
  }
  
  void calcResidual() {
    double infNorm = 0.0;
    
    srand48(thisIndex);
    for (int i = 0; i < n; i++) {
      out[i][0] = out[i][0]/(N*N) - drand48();
      out[i][1] = out[i][1]/(N*N) - drand48();
      
      double mag = sqrt(pow(out[i][0], 2) + pow(out[i][1], 2));
      if (mag > infNorm) infNorm = mag;
    }
    
    double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));
    
    CkCallback cb(CkReductionTarget(main_to_fft, printResidual), m_main);
    contribute(sizeof(double), &r, CkReduction::max_double, cb);
  }
};

#include "fft.def.h"

GCMP_A(fft);
  G_PROPERTY(uint32_t, numChares);
  G_CHARM_APROVIDE(fft_to_main, to_main);
  G_CHARM_USE2(main_to_fft, from_main);
  G_CHARM_AUSE2(transpose, from_transpose)
GEND
