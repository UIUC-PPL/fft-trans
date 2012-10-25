#include <fftw3.h>
#include "TopoManager.h"
#include "NDMeshStreamer.h"

typedef CProxy_GroupChunkMeshStreamer<fftw_complex> streamer_t;

#include "fft1d.decl.h"
PUPbytes(fftw_complex);

#define BUFSIZE 8192 //tunable parameter per machine

#include "data.cc"
#include "fft.cc"

struct Main : public CBase_Main {
  Main_SDAG_CODE

  double start;
  uint64_t N;
  CProxy_fft fftProxy;
  CProxy_fftData data;
  streamer_t streamer;

  Main(CkArgMsg* m) {
    N = atol(m->argv[1]);
    delete m;

    if (N % CkNumPes() != 0)
      CkAbort("CkNumPes() not a factor of N\n");

    // Construct an array of fft chares to do the calculation
    data = CProxy_fftData::ckNew(N);
    init(FFTW_FORWARD, CkCallback(CkReductionTarget(Main,startTimer), thisProxy));
  }

  void init(int sign, CkCallback cb) {
    TopoManager tmgr; // get dimensions for software routing
    int dims[4] = {tmgr.getDimNZ(), tmgr.getDimNY(), tmgr.getDimNX(), tmgr.getDimNT()};
    CkPrintf("Running on NX %d NY %d NZ %d NT %d\n", dims[0], dims[1], dims[2], dims[3]);

    fftProxy = CProxy_fft::ckNew(N, data, sign, cb);
    streamer = streamer_t::ckNew(BUFSIZE, 4, dims, fftProxy);
  }

  void startTimer() {
    start = CkWallTimer();
    // Broadcast the 'go' signal to the fft chare array
    fftProxy.doFFT(CkCallback(CkReductionTarget(Main,stopTimer), thisProxy), streamer);
  }

  void stopTimer() {
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("cores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
             CkNumPes(), N*N, time, gflops);
    validate();
  }
};

#include "fft1d.def.h"
