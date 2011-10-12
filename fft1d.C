/** \file fft1d.C
 *  Author: Abhinav S Bhatele
 *  Date Created: July 16th, 2010
 */

#include "fft1d.decl.h"
#include <math.h>
#include <fftw3.h>

#define TWOPI 6.283185307179586

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ int N;

class fftMsg : public CMessage_fftMsg {
  public:
    int size;
    int source;
    double *data;

    fftMsg() {}

    fftMsg(int _size) : size(_size) {
    }
};

class Main : public CBase_Main {
  public:
    CProxy_fft fftProxy;
    int iteration;

    Main(CkArgMsg* m) {
      numChares = atoi(m->argv[1]);
      N = atoi(m->argv[2]);
      delete m;
      CkPrintf("Each chare will have %d values\n", N*N/numChares);

      CkPrintf("Starting 1D FFT computation ...\n\n");
      fftProxy = CProxy_fft::ckNew(numChares);
      fftProxy.doFFT();
    }

    void done() {
      CkPrintf("All done ...\n");
      CkExit();
    }
};

class fft : public CBase_fft {
  fft_SDAG_CODE

  public:
    int iteration;
    fftw_complex* in; //input data
    fftw_complex* out; //output result
    int n;
    fftw_plan p1;

    fft() {
      __sdag_init();
      iteration = 0;

      n = N*N/(numChares);

      //Initialize data
      //in = new double[n];
      //out = new double[n];
      in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      for (int i=0; i<n; i++) {
        in[i][0] = thisIndex*n+i;
        in[i][1] = 0.0;
        printf("init: [%d].%d = %f\n",thisIndex,i,in[i][0]);
      }
      
      //p1 = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fft(CkMigrateMessage* m) {}
    ~fft() {}

    void sendTranspose()
    {
      //CkPrintf("[%d] sending an array to [%d]\n", thisIndex.x, thisIndex.y, thisIndex.y, thisIndex.x);
      //thisProxy(thisIndex.y,thisIndex.x).getTranspose(real);
      fftMsg **msgs = new fftMsg*[numChares];
      for(int i=0; i<numChares; i++) {
        msgs[i] = new (n/N*2) fftMsg(n/N*2);
        msgs[i]->size = n/N;
        msgs[i]->source = thisIndex;
      }

      int l;
      for(int k=0; k<numChares; k++) {
        l = 0;
        for(int j=0; j<N/numChares; j++) {
          for(int i=0; i<N/numChares; i++) {
            msgs[k]->data[l++] = in[k*N/numChares+(j*N+i)][0];
            msgs[k]->data[l++] = in[k*N/numChares+(j*N+i)][1];
            CkPrintf("[%d].%d %f\n",thisIndex,l-2,in[k*N/numChares+(j*N+i)][0]);
          }
        }
        thisProxy[k].getTranspose(msgs[k]);
      }

      //one message per element
      /*
         for(int j=0; j<N; j++)
         for(int i=j; i<n; i+=N)
         msgs[j]->data[0] = real->data[i];
      //CkPrintf("[%d] packed [%d]=%f for %d\n",thisIndex, i,real->data[i],j);

       */
    }

    void getTranspose(fftMsg *m)
    {
      int k = m->source;
      int l = 0;
      for(int j=0; j<N/numChares; j++)
        for(int i=0; i<N/numChares; i++) {
          out[k*N/numChares+(i*N+j)][0] = m->data[l++];
          out[k*N/numChares+(i*N+j)][1] = m->data[l++];
          CkPrintf("[%d] real[%d] = %f\n",thisIndex,k*N/numChares+(i*N+j),m->data[l-2]);
        }
    }

    void compute()
    {
      fftw_execute(p1); /* repeat as needed */                         
      fftw_destroy_plan(p1);
      CkPrintf("[%d] Computing...\n", thisIndex);
    }

};

#include "fft1d.def.h"
