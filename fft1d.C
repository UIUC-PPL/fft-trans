#include <fftw3.h>
#include <limits>

#define N2 100
#define NCHARE 2
#define BUFSIZE N2*N2/NCHARE/NCHARE

struct fftBuf {
  int source;
  fftw_complex data[BUFSIZE];
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

    mainProxy = thisProxy;
    TopoManager tmgr; // get dimensions for software routing
    int dims[3] = {tmgr.getDimNX()*tmgr.getDimNT(), tmgr.getDimNY(), tmgr.getDimNZ()};
    CkPrintf("Running on NX %d NY %d NZ %d\n", dims[0], dims[1], dims[2]);

    if (N % numChares != 0)
      CkAbort("numChares not a factor of N\n");

    // Construct an array of fft chares to do the calculation
    fftProxy = CProxy_fft::ckNew();
    aggregator = CProxy_GroupMeshStreamer<fftBuf>::ckNew(3, dims, fftProxy, numChares);
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

  int iteration, count, thisIndex;
  uint64_t n;
  fftw_plan p1;
  fftBuf **msgs;
  fftw_complex *in, *out;
  bool validating;

  fft() {
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

    msgs = new fftBuf*[numChares];
    for(int i = 0; i<numChares; i++) msgs[i] = new fftBuf;

    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(Main,FFTReady), mainProxy));
  }

  void sendTranspose(fftw_complex *src) {
    for(int i = 0; i < numChares; i++) {
      msgs[i]->source = thisIndex;
      for(int j = 0, l = 0; j < N/numChares; j++)
        memcpy(msgs[i]->data[(l++)*N/numChares], src[i*N/numChares+j*N], sizeof(fftw_complex)*N/numChares);
    }
    ((GroupMeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->init(1, CkCallback(CkIndex_fft::sendData(), thisProxy), CkCallback(CkIndex_fft::transposeDone(), thisProxy), std::numeric_limits<int>::min(), false);
  }

  void sendData() {
   for(int i = 0; i < numChares; i++)
      ((GroupMeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->insertData(*msgs[i], i);
    ((GroupMeshStreamer<fftBuf> *)CkLocalBranch(aggregator))->done();
  }

  //perform the local transpose
  void process(const fftBuf &m) {
    int k = m.source;
    for(int j = 0, l = 0; j < N/numChares; j++)
      for(int i = 0; i < N/numChares; i++) {
        out[k*N/numChares+(i*N+j)][0] = m.data[l][0];
        out[k*N/numChares+(i*N+j)][1] = m.data[l++][1];
      }
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
