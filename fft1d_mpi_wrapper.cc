#include "fftw3.h"
#include "TopoManager.h"
#include "NDMeshStreamer.h"

#include "fftData.decl.h"
#include "fft1d_mpi_wrapper.decl.h"
#include "mpi-interoperate.h"
PUPbytes(fftw_complex);

//tunable parameter per machine
#define BUFSIZE 131072

fftw_complex *globalDataIn;
fftw_complex *globalDataOut;

struct fftData : public CBase_fftData {
  fftData(CkCallback cb) { contribute(cb); }
  fftw_complex* getIn() { return globalDataIn; }
  fftw_complex* getOut() { return globalDataOut; }
};

#include "fft.cc"
struct Controller : public CBase_Controller {
  double start;
  uint64_t N;
  int sign;
  CProxy_fft fftProxy;
  CProxy_fftData data;
  streamer_t streamer;

  //controller chare constructor creates data group
  Controller(uint64_t _N, int _sign) : N(_N), sign(_sign) {
    data = CProxy_fftData::ckNew(CkCallback(CkReductionTarget(Controller,init), thisProxy));
  }

  //once data group are created, create fft chare-array and streamer
  void init() {
    TopoManager tmgr; // get dimensions for software routing
    int dims[4] = {tmgr.getDimNZ(), tmgr.getDimNY(), tmgr.getDimNX(), tmgr.getDimNT()};
    fftProxy = CProxy_fft::ckNew(N, data, sign, CkCallback(CkReductionTarget(Controller,startTimer),thisProxy));
    streamer = streamer_t::ckNew(BUFSIZE, 4, dims, fftProxy);
  }

  //start the actual fft
  void startTimer() {
    start = CkWallTimer();
    fftProxy.doFFT(CkCallback(CkReductionTarget(Controller,stopTimer), thisProxy), streamer);
  }

  //end of fft
  void stopTimer() {
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("cores: %d\nsize: %ld\ntime: %f sec\nrate: %f GFlop/s\n",
            CkNumPes(), N*N, time, gflops);
    CkExit();
  }
};


//interface function to transfer control to Charm
void fft1d(fftw_complex *_dataIn, fftw_complex *_dataOut, int _N, int sign)
{
  globalDataIn = _dataIn;
  globalDataOut = _dataOut;
  MPI_Barrier(MPI_COMM_WORLD);
  //zero processor creates the controller chare
  if(CkMyPe() == 0) {
    CkPrintf("FFT-1D using Charm++\n");
    CProxy_Controller main = CProxy_Controller::ckNew(_N, sign);
  }
  StartCharmScheduler();
}

#include "fft1d_mpi_wrapper.def.h"
#include "fftData.def.h"
