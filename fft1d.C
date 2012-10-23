#include <fftw3.h>
#include <limits>
#include "charm++.h"

#define N2 100
#define NCHARE 2
#define BUFSIZE N2/NCHARE*N2/NCHARE
#define TOTALBUFSIZE 16384

struct fftBuf {
  int iter;
  int source;
  fftw_complex data[BUFSIZE];

  void pup(PUP::er &p) {
    p | iter;
    p | source;
    PUParray(p,(double*)data,BUFSIZE*2);
  }
};

#include "fileio.h"
#include "TopoManager.h"
#include "NDMeshStreamer.h"
#include "fft1d.decl.h"

#define TWOPI 6.283185307179586

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ uint64_t N;
/*readonly*/ CProxy_GroupMeshStreamer<fftBuf> aggregator;

struct Main : public CBase_Main {
  double start;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
    numChares = CkNumPes();
    N = N2;
    delete m;

    TopoManager tmgr; // get dimensions for software routing
    int dims[4] = {tmgr.getDimNZ(), tmgr.getDimNY(), tmgr.getDimNX(), tmgr.getDimNT()};
    CkPrintf("Running on NX %d NY %d NZ %d NT %d\n", dims[0], dims[1], dims[2], dims[3]);

    mainProxy = thisProxy;

    if (N % numChares != 0)
      CkAbort("numChares not a factor of N\n");

    // Construct an array of fft chares to do the calculation
    fftProxy = CProxy_fft::ckNew();

    CkPrintf("BUFSIZE = %d KB\nTOTALBUFSIZE = %d KB\n", BUFSIZE*16/1024, TOTALBUFSIZE);
    int NUM_MESSAGES_BUF = TOTALBUFSIZE/4/(BUFSIZE*16/1024);
    CkPrintf("NUM_MESSAGES_BUF = %d\nTotal Buf Size = %d KB\n", NUM_MESSAGES_BUF, NUM_MESSAGES_BUF*BUFSIZE*16/1024*4);
    aggregator = CProxy_GroupMeshStreamer<fftBuf>::ckNew(NUM_MESSAGES_BUF, 4, dims, fftProxy);
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

struct fft : public MeshStreamerGroupClient<fftBuf> {
  fft_SDAG_CODE

  int iteration, count;
  uint64_t n;
  fftw_plan p1;
  fftBuf *msg;
  fftw_complex *in, *out;
  bool validating;

  fft() {
    validating = false;

    n = N*N/numChares;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    int length[] = {N};
    p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
                            out, length, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);

    srand48(CkMyPe());
    for(int i = 0; i < n; i++) in[i] = {drand48(), drand48()};

    msg = new fftBuf;
    msg->source = CkMyPe();

    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(Main,FFTReady), mainProxy));
  }

  void initStreamer() {
    ((GroupMeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->init(1, CkCallback(CkIndex_fft::startTranspose(), thisProxy), CkCallback(CkIndex_fft::doneStreaming(), thisProxy), std::numeric_limits<int>::min(), false);
  }

  void sendTranspose(fftw_complex *src_buf) {
    for(int i = 0; i < numChares; i++) {
      for(int j = 0, l = 0; j < N/numChares; j++)
        memcpy(msg->data[(l++)*N/numChares], src_buf[i*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);

      ((GroupMeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->insertData(*msg, i);
    }
    ((GroupMeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->done();
  }

  void applyTranspose(fftBuf &m) {
    for(int j = 0, l = 0; j < N/numChares; j++)
      for(int i = 0; i < N/numChares; i++)
        out[m.source*N/numChares+(i*N+j)] = {m.data[l][0], m.data[l++][1]};
  }

  void twiddle(double sign) {
    double a, c, s, re, im;

    for(int i = 0; i < N/numChares; i++)
      for(int j = 0; j < N; j++) {
        a = sign * (TWOPI*(i+CkMyPe()*N/numChares)*j)/(N*N);
        c = cos(a);
        s = sin(a);

        int idx = i*N+j;

        re = c*out[idx][0] - s*out[idx][1];
        im = s*out[idx][0] + c*out[idx][1];
        out[idx] = {re, im};
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

    srand48(CkMyPe());
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
