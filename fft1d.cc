#include <fftw3_essl.h>
#include "TopoManager.h"
#include "NDMeshStreamer.h"

typedef CProxy_GroupChunkMeshStreamer<fftw_complex> streamer_t;

#include "fft1d.decl.h"
#include "fftData.decl.h"
PUPbytes(fftw_complex);

#include "data.cc"
#include "fft.cc"

struct Main : public CBase_Main {
  Main_SDAG_CODE

  double start;
  uint64_t N;
  CProxy_fft fftProxy;
  CProxy_fftData data;
  streamer_t streamer;
  int NUM_BUF;

  Main(CkArgMsg* m) {
    N = atol(m->argv[1]);
    NUM_BUF = atoi(m->argv[2]);
    delete m;

    if (N % CkNumPes() != 0)
      CkAbort("CkNumPes() not a factor of N\n");

    startBenchmark();
  }

  void init(int sign, CkCallback cb) {
    TopoManager tmgr; // get dimensions for software routing
    int dims[4] = {tmgr.getDimNZ(), tmgr.getDimNY(), tmgr.getDimNX(), tmgr.getDimNT()};
    CkPrintf("Running on NX %d NY %d NZ %d NT %d\n", dims[0], dims[1], dims[2], dims[3]);

    fftProxy = CProxy_fft::ckNew(N, data, sign, cb);
    streamer = streamer_t::ckNew(4, dims, fftProxy, NUM_BUF);
  }
};

#include "fft1d.def.h"
#include "fftData.def.h"
