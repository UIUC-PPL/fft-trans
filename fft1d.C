#include <fftw3.h>
#include <limits>

#define N2 100
#define NCHARE 2
#define BUFSIZE N2*N2/NCHARE/NCHARE

struct fftBuf {
  int iter;
  int source;
  fftw_complex data[BUFSIZE];
};

#include "fileio.h"
#include "TopoManager.h"
#include "MeshStreamer.h"
#include "fft1d.decl.h"

#define TWOPI 6.283185307179586

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ uint64_t N;
/*readonly*/ CProxy_MeshStreamer<fftBuf> aggregator;

struct fftMsg : public CMessage_fftMsg {
  int source;
  fftw_complex *data;
};

struct Main : public CBase_Main {
  double start;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
    numChares = CkNumPes();
    N = N2;
    delete m;

    TopoManager tmgr;
    //use this if you do not want to differentiate based on core ID's
    int NUM_ROWS = tmgr.getDimNX()*tmgr.getDimNT();
    int NUM_COLUMNS = tmgr.getDimNY();
    int NUM_PLANES = tmgr.getDimNZ();
    int NUM_MESSAGES_BUFFERED = numChares;
    CkPrintf("Running on NX %d NY %d NZ %d\n",NUM_ROWS,NUM_COLUMNS,NUM_PLANES);

    mainProxy = thisProxy;

    if (N % numChares != 0)
      CkAbort("numChares not a factor of N\n");

    // Construct an array of fft chares to do the calculation
    fftProxy = CProxy_fft::ckNew();
    aggregator = CProxy_MeshStreamer<fftBuf>::ckNew(NUM_MESSAGES_BUFFERED, NUM_ROWS, NUM_COLUMNS, NUM_PLANES, fftProxy);
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

  void startFlush() {
    //CkPrintf("Starting flush...\n");
    CkStartQD(CkCallback(CkIndex_Main::doFlush(), mainProxy));
  }

  void doFlush() {
    //CkPrintf("Doing flush\n");
    aggregator.flushDirect();
  }

  void printResidual(double r) {
    CkPrintf("residual = %g\n", r);
    CkExit();
  }
};

struct fft : public MeshStreamerClient<fftBuf> {
  fft_SDAG_CODE

  int iteration, count;
  uint64_t n;
  fftw_plan p1;
  fftBuf *msg;
  fftw_complex *in, *out;
  bool validating;
  int thisIndex;

  fft() {
    __sdag_init();

    thisIndex = CkMyPe();
    validating = false;

    n = N*N/numChares;

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

    msg = new fftBuf;

    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(Main,FFTReady), mainProxy));
  }

  void sendTranspose(fftw_complex *src_buf) {
    // All-to-all transpose by constructing and sending
    // point-to-point messages to each chare in the array.
    for(int i = thisIndex; i < thisIndex+numChares; i++) {
      //  Stagger communication order to avoid hotspots and the
      //  associated contention.
      int k = i % numChares;
      for(int j = 0, l = 0; j < N/numChares; j++)
        memcpy(msg->data[(l++)*N/numChares], src_buf[k*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);

      // Tag each message with the iteration in which it was
      // generated, to prevent mis-matched messages from chares that
      // got all of their input quickly and moved to the next step.
      // Runtime system takes ownership of messages once they're sent
      msg->iter = iteration;
      msg->source = thisIndex;
      ((MeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->insertData(msg, k);
    }
  }

  void receiveCombinedData(MeshStreamerMessage<fftBuf> *msg) {
    for(int i = 0; i < msg->numDataItems; i++) {
      fftBuf *m1 = &((msg->data)[i]);
      fftMsg *m2 = new fftMsg;
      m2->source = m1->source;
      m2->data = m1->data;
      CkSetRefNum(m2, m1->iter);
      processData(m2);
    }
  }

  void process(fftBuf* m) {
    /*
    fftMsg *msg = new fftMsg;
    msg->source = m->source;
    memcpy(msg->data, m->data, BUFSIZE*sizeof(fftw_complex));
    //CkPrintf("%d process\n",thisIndex);
    CkSetRefNum(msg, m->iter);
    processData(msg);
    */
  }

  void applyTranspose(fftMsg *m) {
    int k = m->source;
    for(int j = 0, l = 0; j < N/numChares; j++)
      for(int i = 0; i < N/numChares; i++) {
        out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
        out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
      }

    // Save just-received messages to reuse for later sends, to
    // avoid reallocation
    delete m;
    //delete msgs[k];
    //msgs[k] = m;
    //msgs[k]->source = thisIndex;
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

#include "fft1d.def.h"
