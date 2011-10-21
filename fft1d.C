/** \file fft1d.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 16th, 2010
 */

#include "fft1d.decl.h"
#include <math.h>
#include <fftw3.h>
#include "verify.h"
#include "register.h"

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
	       N, numChares, time, 5*N*N*log2(N*N)/(time*1000000000));
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
    fftMsg **msgs;

    double copytime, computetime, sendtime, twiddletime;

    fft() {
      __sdag_init();

      n = N*N/(numChares);

      copytime = computetime = sendtime = 0.0;

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
          out, inembed,
          istride, idist,
          out, onembed,
          ostride, odist,
          FFTW_FORWARD, FFTW_ESTIMATE);

      //TODO: numbers need to be generated independantly?
      srand48(thisIndex);
      for (int i=0; i<n; i++) {
        in[i][0] = drand48();
        in[i][1] = drand48();
      }

      msgs = new fftMsg*[numChares];
      for(int i=0; i<numChares; i++) {
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
      //CkPrintf("[%d] sending an array to [%d]\n", thisIndex.x, thisIndex.y, thisIndex.y, thisIndex.x);
      //thisProxy(thisIndex.y,thisIndex.x).getTranspose(real);
      fftw_complex *buf = (iteration == 0) ? in : out;

      double this_copytime = CkWallTimer();

      for(int k=0; k<numChares; k++) {
        int l = 0;
        envelope *env = UsrToEnv(msgs[k]);
        _SET_USED(env, 0);
        if (env->isPacked()) {
          //CkPrintf("[%d] isPacked\n", thisIndex);
          unsigned char msgidx = env->getMsgIdx();
          if(_msgTable[msgidx]->unpack) {
            msgs[k] = (fftMsg*)_msgTable[msgidx]->unpack(msgs[k]);
            UsrToEnv(msgs[k])->setPacked(0);
          }
        }

        for(int j=0; j<N/numChares; j++)
          memcpy(msgs[k]->data[(l++)*N/numChares], buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);
        //thisProxy[k].getTranspose(msgs[k+offset]);
      }

      copytime += CkWallTimer() - this_copytime;

      double this_sendtime = CkWallTimer();

      for(int k=0; k<numChares; k++) {
        CmiReference(UsrToEnv(msgs[k]));
        thisProxy[k].getTranspose(msgs[k]);
      }
    
      sendtime += CkWallTimer() - this_sendtime;
    }

    void applyTranspose(fftMsg *m)
    {
      int k = m->source;
      int l = 0;
      for(int j=0; j<N/numChares; j++)
        for(int i=0; i<N/numChares; i++) {
          out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
          out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
          //CkPrintf("[%d] real[%d] = %f\n",thisIndex,k*N/numChares+(i*N+j),m->data[l-2]);
        }
    }

    void compute(bool doTwiddle)
    {
      double this_computetime = CkWallTimer();

      fftw_execute(p1);

      //for(int i=0; i<n; i++)
      //CkPrintf("[%d] in[%d] = %f -> out[%d] = %f\n",thisIndex, i, in[i][0], i, out[i][0]);
      //CkPrintf("[%d] in[%d] = %f + %fi\n",thisIndex, i, in[i][0], in[i][1]);

      //fftw_destroy_plan(p1);
      //CkPrintf("[%d] Computing...\n", thisIndex);

      if(doTwiddle) {
        twiddletime = CkWallTimer();
        twiddle();
        twiddletime = CkWallTimer() - twiddletime;
      }

      computetime += CkWallTimer() - this_computetime;
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

          re = c*out[idx][0] - s*out[idx][1];
          im = s*out[idx][0] + c*out[idx][1];
          out[idx][0] = re;
          out[idx][1] = im;
        }
    }
};

#include "fft1d.def.h"
