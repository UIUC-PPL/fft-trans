#include <fftw3.h>
#include "TopoManager.h"
#include "NDMeshStreamer.h"

typedef CProxy_GroupChunkMeshStreamer<fftw_complex> streamer_t;

#include "fft1d_mpi_wrapper.decl.h"
#include "fftData.decl.h"
#include "mpi-interoperate.h"
PUPbytes(fftw_complex);

#define BUFSIZE 8192 //tunable parameter per machine


fftw_complex *globalDataIn;
fftw_complex *globalDataOut;

struct fftData : public CBase_fftData {
  fftData(CkCallback cb) { contribute(cb); }
  fftw_complex* getIn() { return globalDataIn; }
  fftw_complex* getOut() { return globalDataOut; }
  void swap(CkCallback cb) { fftw_complex *tmp; tmp = globalDataIn; globalDataIn = globalDataOut; globalDataOut = tmp; contribute(cb); }
};

#include "fft.cc"
struct Main : public CBase_Main {
  Main_SDAG_CODE

  double start;
  uint64_t N;
  CProxy_fft fftProxy;
  CProxy_fftData data;
  streamer_t streamer;

  Main(uint64_t _N) {
    N = _N;
    if (N % CkNumPes() != 0)
      CkAbort("CkNumPes() not a factor of N\n");

    startBenchmark();
  }

  void init(int sign, CkCallback cb) {
    TopoManager tmgr; // get dimensions for software routing
    int dims[4] = {tmgr.getDimNZ(), tmgr.getDimNY(), tmgr.getDimNX(), tmgr.getDimNT()};
    CkPrintf("Running on NX %d NY %d NZ %d NT %d\n", dims[0], dims[1], dims[2], dims[3]);

    fftProxy = CProxy_fft::ckNew(N, data, sign, cb);
    streamer = streamer_t::ckNew(BUFSIZE, 4, dims, fftProxy);
  }
};

void fft1d(fftw_complex *_dataIn, fftw_complex *_dataOut, int _N)
{
  globalDataIn = _dataIn;
  globalDataOut = _dataOut;
  MPI_Barrier(MPI_COMM_WORLD);
  if(CkMyPe() == 0) {
    CkPrintf("FFT-1D using Charm++\n");
    CProxy_Main main = CProxy_Main::ckNew(_N);
  }
  CsdScheduler(-1);
}



#include "fft1d_mpi_wrapper.def.h"
#include "fftData.def.h"
