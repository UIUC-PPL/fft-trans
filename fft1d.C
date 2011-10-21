/** \file fft1d.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 16th, 2010
 */

#include "fft1d.decl.h"
#include <fftw3.h>
#include "verify.h"

#define TWOPI 6.283185307179586

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ int N;

struct fftMsg : public CMessage_fftMsg {
    int source;
    fftw_complex *data;
};

void printMat(fftw_complex *data, int x, int y, const char* msg, int idx) {
  CkPrintf("[%d] %s\n",idx,msg);
  for(int i=0; i<x; i++) {
    for(int j=0; j<y; j++) {
     CkPrintf(" %5.0f", data[i*y+j][0]);
    }
    CkPrintf("\n");
  }
}

struct Main : public CBase_Main {
  double start;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
      numChares = atoi(m->argv[1]);
      N = atoi(m->argv[2]);
      delete m;

      mainProxy = thisProxy;

      fftProxy = CProxy_fft::ckNew(numChares);
    }

  void startTiming() {
      start = CkWallTimer();
      fftProxy.doFFT();
  }

    void done() {
      double time = CkWallTimer() - start;
      CkPrintf("FFT on %d elements with %d chares took %f time at %f GFlop/s\n",
                N, numChares, time, 5*N*N*log2(N*N)/(time*1000000000));
      fftProxy.writeResults(0);
    }
};

struct fft : public CBase_fft {
  fft_SDAG_CODE

    int iteration, count;
    int n;
    fftw_plan p1;
    fftMsg **msgs;
    fftw_complex *in, *out;

    fft() {
      __sdag_init();

      n = N*N/numChares;

      in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

      int rank = 1; /* not 2: we are computing 1d transforms */
      int length[] = {N}; /* 1d transforms of length 10 */
      int howmany = N/numChares;
      int idist = N;
      int odist = N;
      int istride = 1;
      int ostride = 1;
      int *inembed = length, *onembed = length;

      p1 = fftw_plan_many_dft(rank, length, howmany,
              out, inembed,
              istride, idist,
              out, onembed,
              ostride, odist,
              FFTW_FORWARD, FFTW_ESTIMATE);

      srand48(thisIndex);
      for(int i=0; i<n; i++) {
        in[i][0] = drand48();
        in[i][1] = drand48();
      }

      //printMat(buf, N/numChares, N, "Initialized", thisIndex);

      msgs = new fftMsg*[numChares*3];
      for(int i=0; i<numChares*3; i++) {
        msgs[i] = new (n/numChares) fftMsg;
        msgs[i]->source = thisIndex;
      }

      contribute(CkCallback(CkIndex_Main::startTiming(), mainProxy));
    }

    fft(CkMigrateMessage* m) {}
    ~fft() {}

    void sendTranspose()
    {
      if(thisIndex == 0)
        CkPrintf("TRANSPOSING\n");

      fftw_complex *buf = (iteration == 0) ? in : out;

      int offset = iteration*numChares;

      for(int k=0; k<numChares; k++) {

        int l = 0;
        for(int j=0; j<N/numChares; j++)
          memcpy(msgs[k+offset]->data[(l++)*N/numChares], buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);
        //thisProxy[k].getTranspose(msgs[k+offset]);
      }

      for(int k=0; k<numChares; k++) {
        CkSetRefNum(msgs[k+offset], iteration);
        thisProxy[k].getTranspose(msgs[k+offset]);
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
    }

    void compute(bool doTwiddle)
    {
      //printMat(buf, N/numChares, N, "AtCompute", thisIndex);

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

#include "fft1d.def.h"
