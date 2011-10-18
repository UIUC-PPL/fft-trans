/** \file fft1d.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 16th, 2010
 */

#include "fft1d.decl.h"
#include <math.h>
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

struct Main : public CBase_Main {
  double start;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
      numChares = atoi(m->argv[1]);
      N = atoi(m->argv[2]);
      delete m;

      if(N%numChares !=0)
        CkAbort("numChares not a multiple of N\n");

      fftProxy = CProxy_fft::ckNew(numChares);
    }

  void startTiming() {
      start = CkWallTimer();
      fftProxy.doFFT();
  }

    void done() {
      double time = CkWallTimer() - start;
      CkPrintf("FFT on %d elements with %d chares took %f time at %f GFlop/s\n",
	       N, numChares, time, 5*N*log2(N)/(time*1000000000));
      fftProxy.writeResults(0);
    }
};

struct fft : public CBase_fft {
  fft_SDAG_CODE

    int iteration, count;
    fftw_complex* in; //input data
    fftw_complex* out; //output result
    int n;
    fftw_plan p1;

    fft() {
      __sdag_init();

      n = N*N/(numChares);

      //Initialize data
      //in = new double[n];
      //out = new double[n];
      in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

      int rank = 1; /* not 2: we are computing 1d transforms */
      int length[] = {N}; /* 1d transforms of length 10 */
      int howmany = N/numChares;
      int idist = N;
      int odist = N;
      int istride = 1;
      int ostride = 1; /* distance between two elements in
                                    the same column */
      int *inembed = length, *onembed = length;

      p1 = fftw_plan_many_dft(rank, length, howmany,
          in, inembed,
          istride, idist,
          in, onembed,
          ostride, odist,
          FFTW_FORWARD, FFTW_ESTIMATE);

      //TODO: numbers need to be generated independantly?
      srand48(thisIndex);
      for (int i=0; i<n; i++) {
        in[i][0] = drand48();
        in[i][1] = drand48();
      }

      contribute(CkCallback(CkIndex_Main::startTiming(), mainProxy));
    }

    fft(CkMigrateMessage* m) {}
    ~fft() {}

    void sendTranspose()
    {
      if(thisIndex == 0)
        CkPrintf("TRANSPOSING\n");
      //CkPrintf("[%d] sending an array to [%d]\n", thisIndex.x, thisIndex.y, thisIndex.y, thisIndex.x);
      //thisProxy(thisIndex.y,thisIndex.x).getTranspose(real);
      fftMsg **msgs = new fftMsg*[numChares];
      for(int i=0; i<numChares; i++) {
        msgs[i] = new (n/numChares) fftMsg;
        msgs[i]->source = thisIndex;
      }

      for(int k=0; k<numChares; k++) {
        int l = 0;
        for(int j=0; j<N/numChares; j++) {
          for(int i=0; i<N/numChares; i++) {
            msgs[k]->data[l][0] = in[k*N/numChares+(j*N+i)][0];
            msgs[k]->data[l++][1] = in[k*N/numChares+(j*N+i)][1];
            //CkPrintf("[%d].%d %f\n",thisIndex,l-2,in[k*N/numChares+(j*N+i)][0]);
          }
        }
        thisProxy[k].getTranspose(msgs[k]);
      }
    }

    void applyTranspose(fftMsg *m)
    {
      int k = m->source;
      int l = 0;
      for(int j=0; j<N/numChares; j++)
        for(int i=0; i<N/numChares; i++) {
          in[k*N/numChares+(i*N+j)][0] = m->data[l][0];
          in[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
          //CkPrintf("[%d] real[%d] = %f\n",thisIndex,k*N/numChares+(i*N+j),m->data[l-2]);
        }
    }

    void compute(bool doTwiddle)
    {
      fftw_execute(p1);

      //for(int i=0; i<n; i++)
      //CkPrintf("[%d] in[%d] = %f -> out[%d] = %f\n",thisIndex, i, in[i][0], i, out[i][0]);
      //CkPrintf("[%d] in[%d] = %f + %fi\n",thisIndex, i, in[i][0], in[i][1]);

      //fftw_destroy_plan(p1);
      //CkPrintf("[%d] Computing...\n", thisIndex);

      if(doTwiddle)
        twiddle();
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

          //CkPrintf("[%d] Twiddle for [%d,%d]: tw_re = %f tw_im = %f\n",thisIndex,(i+k*N/numChares), j, c, s);

          re = c*in[idx][0] - s*in[idx][1];
          im = s*in[idx][0] + c*in[idx][1];
          in[idx][0] = re;
          in[idx][1] = im;
        }
    }
};

#include "fft1d.def.h"
