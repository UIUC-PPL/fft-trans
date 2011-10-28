#include "fft1d.decl.h"
#include <fftw3.h>
#include <limits>
#include "fileio.h"

#define TWOPI 6.283185307179586

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ uint64_t N, n;
/*readonly*/ CProxy_transpose transposer;

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

    n = N*N/numChares;
    mainProxy = thisProxy;

    if (N % numChares != 0)
      CkAbort("numChares not a factor of N\n");

    // Construct an array of fft chares to do the calculation
    fftProxy = CProxy_fft::ckNew(numChares);

    CkArrayOptions opts(numChares);
    opts.bindTo(fftProxy);
    transposer = CProxy_p2p_transpose::ckNew(opts);
  }

  void FFTReady() {
    start = CkWallTimer();
    // Broadcast the 'go' signal to the fft chare array
    fftProxy.doFFT();
  }

  void FFTDone() {
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("chares: %d\ncores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
             numChares, CkNumPes(), N*N, time, gflops);

    fftProxy.initValidation();
  }

  void printResidual(double r) {
    CkPrintf("residual = %g\n", r);
    CkExit();
  }
};

struct transpose : public CBase_transpose
{
  transpose() { }
  transpose(CkMigrateMessage *) { }
  virtual void sendTranspose(int, fftw_complex*, fft *fft_obj) { CkAbort("Should not be here"); }
};

struct fft : public CBase_fft {
  fft_SDAG_CODE

  int iteration, count;
  fftw_plan p1;
  fftw_complex *in, *out;
  bool validating;

  fft() {
    __sdag_init();

    validating = false;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    int length[] = {N};
    p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
                            out, length, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);

    srand48(thisIndex);
    for(int i = 0; i < n; i++) {
      in[i][0] = drand48();
      in[i][1] = drand48();
    }

    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(Main,FFTReady), mainProxy));
  }

  void twiddle(double sign) {
    double a, c, s, re, im;

    int k = thisIndex;
    for(int i = 0; i < N/numChares; i++)
      for(int j = 0; j < N; j++) {
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

    contribute(CkCallback(CkReductionTarget(Main,FFTReady), mainProxy));
  }

  void calcResidual() {
    double infNorm = 0.0;

    srand48(thisIndex);
    for(int i = 0; i < n; i++) {
      out[i][0] = out[i][0]/(N*N) - drand48();
      out[i][1] = out[i][1]/(N*N) - drand48();

      double mag = sqrt(pow(out[i][0], 2) + pow(out[i][1], 2));
      if(mag > infNorm) infNorm = mag;
    }

    double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));

    CkCallback cb(CkReductionTarget(Main, printResidual), mainProxy);
    contribute(sizeof(double), &r, CkReduction::max_double, cb);
  }

  fft(CkMigrateMessage* m) {}
  ~fft() {}
};

struct p2p_transpose : public CBase_p2p_transpose
{
  fft *fft_obj;
  fftMsg **msgs;
  int count;

  p2p_transpose() {
    __sdag_init();
    msgs = new fftMsg*[numChares];
    for(int i = 0; i < numChares; i++) {
      msgs[i] = new (n/numChares) fftMsg;
      msgs[i]->source = thisIndex;
    }
  }
  p2p_transpose(CkMigrateMessage *) { }

  p2p_transpose_SDAG_CODE

  void sendTranspose(int iteration, fftw_complex *src_buf, fft *obj) {
    fft_obj = obj;

    // All-to-all transpose by constructing and sending
    // point-to-point messages to each chare in the array.
    for(int i = thisIndex; i < thisIndex+numChares; i++) {
      //  Stagger communication order to avoid hotspots and the
      //  associated contention.
      int k = i % numChares;
      for(int j = 0, l = 0; j < N/numChares; j++)
        memcpy(msgs[k]->data[(l++)*N/numChares], src_buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);

      // Tag each message with the iteration in which it was
      // generated, to prevent mis-matched messages from chares that
      // got all of their input quickly and moved to the next step.
      CkSetRefNum(msgs[k], iteration);
      thisProxy[k].getTranspose(msgs[k]);
      // Runtime system takes ownership of messages once they're sent
      msgs[k] = NULL;
    }

    thisProxy[thisIndex].receive(iteration);
  }

  void applyTranspose(fftMsg *m, fftw_complex *out) {
    int k = m->source;
    for(int j = 0, l = 0; j < N/numChares; j++)
      for(int i = 0; i < N/numChares; i++) {
        out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
        out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
      }

    // Save just-received messages to reuse for later sends, to
    // avoid reallocation
    delete msgs[k];
    msgs[k] = m;
    msgs[k]->source = thisIndex;
  }
};

#include "fft1d.def.h"
