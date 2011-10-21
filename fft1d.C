/** \file fft1d.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 16th, 2010
 */

typedef double fftw_complex[2];

#include "fft1d.decl.h"

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
      fftProxy.writeResults(0);
    }
};

struct fft : public CBase_fft {
  fft_SDAG_CODE

    int iteration, count;
    int n;
    fftMsg **msgs;
    fftw_complex *in, *out;

    fft() {
      __sdag_init();

      n = N*N/numChares;

      in = new fftw_complex[n];
      for(int i=0; i<n; i++) {
        in[i][0] = i+thisIndex*n;
        in[i][1] = -i-thisIndex*n;
      }
      out = new fftw_complex[n];

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

      int n = N*N/numChares;
      if(iteration == 1) {
        for(int i=0; i<n; i++) {
          if(out[i][0] != i+thisIndex*n || out[i][1] != -i-thisIndex*n)
            CkPrintf("ERROR: transpose failed\n");
        }
      }
    }
};

#include "fft1d.def.h"
