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

class fftMsg : public CMessage_fftMsg {
  public:
    int size;
    int source;
    fftw_complex *data;

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

      if(N%numChares !=0)
        CkAbort("numChares not a multiple of N\n");

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
    fftw_plan* plans;
    int count;
    fftw_plan p1;

    int transposeCount;
    int computeCount;

    fft() {
      __sdag_init();
      iteration = 0;

      count = 0;
      transposeCount = 0;
      computeCount = 0;

      n = N*N/(numChares);

      //Initialize data
      //in = new double[n];
      //out = new double[n];
      in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

      fftw_plan* plans = new fftw_plan[N/numChares];
      for(int i=0; i<N/numChares; i++)
        plans[i] = fftw_plan_dft_1d(N, &in[i*N], &in[i*N], FFTW_FORWARD, FFTW_ESTIMATE);

      p1 = fftw_plan_dft_1d(N, &in[0], &in[0], FFTW_FORWARD, FFTW_ESTIMATE);

      for (int i=0; i<n; i++) {
        in[i][0] = thisIndex*n+i;
        in[i][1] = 0.0;
        printf("init: [%d].%d = %f\n",thisIndex,i,in[i][0]);
      }
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
        msgs[i] = new (n/numChares) fftMsg(n/numChares);
        msgs[i]->source = thisIndex;
      }

      int l;
      for(int k=0; k<numChares; k++) {
        l = 0;
        for(int j=0; j<N/numChares; j++) {
          for(int i=0; i<N/numChares; i++) {
            msgs[k]->data[l][0] = in[k*N/numChares+(j*N+i)][0];
            msgs[k]->data[l++][1] = in[k*N/numChares+(j*N+i)][1];
            //CkPrintf("[%d].%d %f\n",thisIndex,l-2,in[k*N/numChares+(j*N+i)][0]);
          }
        }
        thisProxy[k].getTranspose(msgs[k]);
      }

      transposeCount++;
    }

    void getTranspose(fftMsg *m)
    {
      int k = m->source;
      int l = 0;
      for(int j=0; j<N/numChares; j++)
        for(int i=0; i<N/numChares; i++) {
          in[k*N/numChares+(i*N+j)][0] = m->data[l][0];
          in[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
          //CkPrintf("[%d] real[%d] = %f\n",thisIndex,k*N/numChares+(i*N+j),m->data[l-2]);
        }

      count++;

      if(transposeCount == 3){
        contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_Main::done(), mainProxy));
        char filename[80];
        sprintf(filename,"%d-%d.dump%d",numChares,N,thisIndex);
        writeCommFile(n, in, filename);
      }
      else if(count == numChares)
        compute();
    }

    void compute()
    {
      count = 0;
      for(int i=0; i<N/numChares; i++)
        fftw_execute_dft(p1,&in[i*N],&in[i*N]);
      //fftw_execute(plans[0]);

      for(int i=0; i<n; i++)
        //CkPrintf("[%d] in[%d] = %f -> out[%d] = %f\n",thisIndex, i, in[i][0], i, out[i][0]);
        CkPrintf("[%d] in[%d] = %f + %fi\n",thisIndex, i, in[i][0], in[i][1]);

      //fftw_destroy_plan(p1);
      CkPrintf("[%d] Computing...\n", thisIndex);
      if(computeCount == 0)
        twiddle();
      computeCount++;

      if(thisIndex == 0)
        thisProxy.sendTranspose();
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

            CkPrintf("[%d] Twiddle for [%d,%d]: tw_re = %f tw_im = %f\n",thisIndex,(i+k*N/numChares), j, c, s);

            re = c*in[idx][0] - s*in[idx][1];
            im = s*in[idx][0] + c*in[idx][1];
            in[idx][0] = re;
            in[idx][1] = im;
          }
    }
};

#include "fft1d.def.h"
