#include <fftw3.h>
#include "fileio.h"
#include "TopoManager.h"
#include "NDMeshStreamer.h"
#include "fft1d.decl.h"
PUPbytes(fftw_complex);

#define BUFSIZE 8192 //tunable parameter per machine
#define TWOPI 6.283185307179586

/*readonly*/ int numChares;
/*readonly*/ CProxy_GroupChunkMeshStreamer<fftw_complex> aggregator;

struct Main : public CBase_Main {
  double start;
  uint64_t N;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
    numChares = CkNumPes();
    N = atol(m->argv[1]);
    delete m;

    TopoManager tmgr; // get dimensions for software routing
    int dims[4] = {tmgr.getDimNZ(), tmgr.getDimNY(), tmgr.getDimNX(), tmgr.getDimNT()};
    CkPrintf("Running on NX %d NY %d NZ %d NT %d\n", dims[0], dims[1], dims[2], dims[3]);

    if (N % numChares != 0)
      CkAbort("numChares not a factor of N\n");

    // Construct an array of fft chares to do the calculation
    fftProxy = CProxy_fft::ckNew(N, thisProxy);
    aggregator = CProxy_GroupChunkMeshStreamer<fftw_complex>::ckNew(BUFSIZE, 4, dims, fftProxy);
  }

  void FFTReady() {
    start = CkWallTimer();
    // Broadcast the 'go' signal to the fft chare array
    fftProxy.doFFT(CkCallback(CkReductionTarget(Main,FFTDone), thisProxy));
  }

  void FFTDone() {
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("chares: %d\ncores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
             numChares, CkNumPes(), N*N, time, gflops);

    fftProxy.initValidation();
  }
};

#define SET_VALUES(a,b,c)  do { (a)[0] = b; (a)[1] = c; } while (0);

struct fft : public MeshStreamerGroupClient<fftw_complex> {
  fft_SDAG_CODE

  int iteration, count;
  uint64_t n, N;
  fftw_plan p1;
  fftw_complex *in, *out, *buf;
  bool validating;

  fft(uint64_t N, CProxy_Main mainProxy) : validating(false), n(N/CkNumPes()), N(N) {
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n*N);

    int length[] = {(int)N};
    p1 = fftw_plan_many_dft(1, length, n, out, length, 1, N,
                            out, length, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);

    srand48(CkMyPe());
    for(int i = 0; i < n*N; i++) SET_VALUES(in[i], drand48(), drand48());

    buf = new fftw_complex[n*n];

    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(Main,FFTReady), mainProxy));
  }

  void initStreamer() {
    aggregator.ckLocalBranch()->init(1, CkCallback(CkIndex_fft::streamerReady(), thisProxy), CkCallback(CkIndex_fft::doneStreaming(), thisProxy), 0, false);
  }

  void sendTranspose(fftw_complex *src_buf) {
    for(int i = 0; i < numChares; i++) {
      for(int j = 0, l = 0; j < n; j++)
        memcpy(buf[(l++)*n], src_buf[i*n+j*N], sizeof(fftw_complex)*n);

      aggregator.ckLocalBranch()->insertData(buf, n*n, i);
    }
    aggregator.ckLocalBranch()->done();
  }

  void applyTranspose(fftw_complex *data, int numItems, int src) {
    for(int j = 0, l = 0; j < n; j++)
      for(int i = 0; i < n; i++)
        SET_VALUES(out[src*n+(i*N+j)], data[l][0], data[l++][1]);
  }

  void twiddle(double sign) {
    double a, c, s, re, im;

    for(int i = 0; i < n; i++)
      for(int j = 0; j < N; j++) {
        a = sign * (TWOPI*(i+CkMyPe()*n)*j)/(N*N);
        c = cos(a);
        s = sin(a);

        int idx = i*N+j;

        re = c*out[idx][0] - s*out[idx][1];
        im = s*out[idx][0] + c*out[idx][1];
        SET_VALUES(out[idx], re, im);
      }
  }

  void initValidation();
  void calcResidual();
  void printResidual(double residual);
};

#include "verify.cc"
#include "fft1d.def.h"
