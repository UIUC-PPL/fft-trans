#include "fft1d.decl.h"
#include <fftw3.h>
#include <limits>
#include "verify.h"

#define TWOPI 6.283185307179586

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ uint64_t N;

struct fftMsg : public CMessage_fftMsg {
    int source;
    fftw_complex *data;
};

struct Main : public CBase_Main {
  double start;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
      numChares = atoi(m->argv[1]);
      N = atol(m->argv[2]);
      delete m;

      mainProxy = thisProxy;

      if(N%numChares !=0)
        CkAbort("numChares not a multiple of N\n");

      fftProxy = CProxy_fft::ckNew(numChares);
    }

  void startFFT() {
      start = CkWallTimer();
      fftProxy.doFFT();
  }

  void doneFFT() {
      double time = CkWallTimer() - start;
      double gflops = 5*(double)N*N*log2((double)N*N)/(time*1000000000);
      CkPrintf("chares: %d\ncores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
        numChares, CkNumPes(), N*N, time, gflops);

      fftProxy.initValidation();
  }

  void printResidual(CkReductionMsg *m) {
    double *r = (double *)m->getData();
    CkPrintf("residual = %g\n", *r);
    CkExit();
  }
};

struct fft : public CBase_fft {
  fft_SDAG_CODE

    int iteration, count;
    uint64_t n;
    fftw_plan p1;
    fftMsg **msgs;
    fftw_complex *in, *out;
    bool validating;

    fft() {
      __sdag_init();

      validating = false;

      n = N*N/numChares;

      in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

      int length[] = {N};
      p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
              out, length, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);

      srand48(thisIndex);
      for(int i=0; i<n; i++) {
        in[i][0] = drand48();
        in[i][1] = drand48();
      }

      msgs = new fftMsg*[numChares];
      for(int i=0; i<numChares; i++) {
        msgs[i] = new (n/numChares) fftMsg;
        msgs[i]->source = thisIndex;
      }

      contribute(CkCallback(CkIndex_Main::startFFT(), mainProxy));
    }

    fft(CkMigrateMessage* m) {}
    ~fft() {}

    void sendTranspose()
    {
      if(thisIndex == 0)
        CkPrintf("TRANSPOSING\n");

      fftw_complex *buf = (iteration == 0) ? in : out;

      for(int i=thisIndex; i<thisIndex+numChares; i++) {

        int l = 0;
        int k = i % numChares;
        for(int j=0; j<N/numChares; j++)
          memcpy(msgs[k]->data[(l++)*N/numChares], buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);

        CkSetRefNum(msgs[k], iteration);
        thisProxy[k].getTranspose(msgs[k]);
        msgs[k] = NULL;
      }
    }

    void applyTranspose(fftMsg *m)
    {
      int k = m->source;
      int l = 0;
      for(int j=0; j<N/numChares; j++)
        for(int i=0; i<N/numChares; i++) {
          out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
          out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
        }
      delete msgs[k];
      msgs[k] = m;
      msgs[k]->source = thisIndex;
    }

    void compute(bool doTwiddle)
    {
      fftw_execute(p1);
      if(doTwiddle) {
        twiddle();
      }
    }

    void twiddle() {
      double a, c, s, re, im;

      int k = thisIndex;
      for(int i = 0; i<N/numChares; i++)
        for( int j = 0; j<N; j++) {
          a = -(TWOPI*(i+k*N/numChares)*j)/(N*N);
          if(validating) a *= -1;
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

      contribute(CkCallback(CkIndex_Main::startFFT(), mainProxy));
    }

  void calcResidual() {
    double infNorm = 0.0;

    srand48(thisIndex);
    for(int i=0; i<n; i++) {
      out[i][0] = out[i][0]/(N*N) - drand48();
      out[i][1] = out[i][1]/(N*N) - drand48();

      double mag = sqrt(pow(out[i][0],2) + pow(out[i][1], 2));
      if(mag > infNorm) infNorm = mag;
    }

    double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));

    CkCallback cb(CkIndex_Main::printResidual(NULL), mainProxy);
    contribute(sizeof(double), &r, CkReduction::max_double, cb);
  }
};

#include "fft1d.def.h"
