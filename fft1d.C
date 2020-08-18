#include "fft1d.decl.h"
#include "fft1d.h"
#include <limits>
#include "fileio.h"

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int numChares;
/*readonly*/ int N;
/*readonly*/ bool validate;

#ifdef MODE_CUDA
extern void invokeTwiddle(complex_t* out, int N, int numChares, int k,
    double sign, cudaStream_t stream);
#endif

struct fftMsg : public CMessage_fftMsg {
  int source;
  complex_t *data;
};

struct Main : public CBase_Main {
  double start;
  CProxy_fft fftProxy;

  Main(CkArgMsg* m) {
    // Default parameters
    numChares = 4;
    N = 128;
    validate = false;

    int c;
    while ((c = getopt(m->argc, m->argv, "c:n:v")) != -1 ) {
      switch (c) {
        case 'c':
          numChares = atoi(optarg);
          break;
        case 'n':
          N = atoi(optarg);
          break;
        case 'v':
          validate = true;
          break;
        default:
          CkAbort("Invalid parameter");
          break;
      }
    }
    delete m;

    mainProxy = thisProxy;

    if (N % numChares != 0) {
      CkAbort("numChares not a factor of N\n");
    }

    CkPrintf("\n1D FFT\n");
    CkPrintf("\tChares: %d\n\tN: %d, Size: %d\n", numChares, N, N*N);
    CkPrintf("\tValidate: %s\n\n", validate ? "true" : "false");

    // Construct an array of FFT chares to do the calculation
    fftProxy = CProxy_fft::ckNew(numChares);
  }

  void FFTReady() {
    start = CkWallTimer();
    // Broadcast the 'go' signal to the fft chare array
    fftProxy.doFFT();
  }

  void FFTDone() {
    double time = CkWallTimer() - start;
    double gflops = 5 * (double)N*N * log2((double)N*N) / (time * 1000000000);
    CkPrintf("\nTime: %f sec\nRate: %f GFlop/s\n", time, gflops);

    if (validate) {
      CkPrintf("\nPerforming validation\n");
      fftProxy.initValidation();
    } else {
      CkExit();
    }
  }

  void printResidual(double r) {
    CkPrintf("\nResidual = %g\n", r);
    CkExit();
  }
};

struct fft : public CBase_fft {
  fft_SDAG_CODE

  int iteration, count;
  uint64_t n;
#ifdef MODE_CPU
  complex_t *in, *out;
  fftw_plan p1;
#elif defined MODE_CUDA
  complex_t *h_in, *h_out;
  complex_t *d_in, *d_out;
  cufftHandle p1;
  cudaStream_t compute_stream;
  cudaStream_t comm_stream;
  cudaEvent_t compute_event;
  cudaEvent_t comm_event;
#endif
  fftMsg **msgs;
  bool validating;

  fft() {
    __sdag_init();

    validating = false;

    n = N*N/numChares;

#ifdef MODE_CPU
    in = (complex_t*) fftw_malloc(sizeof(complex_t) * n);
    out = (complex_t*) fftw_malloc(sizeof(complex_t) * n);
#elif defined MODE_CUDA
    hapiCheck(cudaMallocHost(&h_in, sizeof(complex_t) * n));
    hapiCheck(cudaMallocHost(&h_out, sizeof(complex_t) * n));
    hapiCheck(cudaMalloc(&d_in, sizeof(complex_t) * n));
    hapiCheck(cudaMalloc(&d_out, sizeof(complex_t) * n));

    hapiCheck(cudaStreamCreateWithPriority(&compute_stream, cudaStreamDefault, 0));
    hapiCheck(cudaStreamCreateWithPriority(&comm_stream, cudaStreamDefault, -1));

    hapiCheck(cudaEventCreateWithFlags(&compute_event, cudaEventDisableTiming));
    hapiCheck(cudaEventCreateWithFlags(&comm_event, cudaEventDisableTiming));
#endif

    int length[] = {N};
#ifdef MODE_CPU
    p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
                            out, length, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);
#elif defined MODE_CUDA
    if (cufftPlanMany(&p1, 1, length, length, 1, N, length, 1, N, CUFFT_C2C,
          N/numChares) != CUFFT_SUCCESS) {
      CkAbort("CUFFT Error: Unable to create plan");
    }
    if (cufftSetStream(p1, compute_stream) != CUFFT_SUCCESS) {
      CkAbort("CUFFT Error: Unable to set CUDA stream");
    }
#endif

    srand48(thisIndex);
#ifdef MODE_CPU
    for (int i = 0; i < n; i++) {
      in[i][0] = drand48();
      in[i][1] = drand48();
    }
#elif defined MODE_CUDA
    // FIXME: Random number generation and assignment happens on the host
    for (int i = 0; i < n; i++) {
      h_in[i].x = static_cast<float>(drand48());
      h_in[i].y = static_cast<float>(drand48());
    }
    hapiCheck(cudaMemcpyAsync(d_in, h_in, sizeof(complex_t) * n,
          cudaMemcpyHostToDevice, comm_stream));
    CkCallback* cb = new CkCallback(CkIndex_fft::initComplete(), thisProxy[thisIndex]);
    hapiAddCallback(comm_stream, cb);
#endif

    msgs = new fftMsg*[numChares];
    for (int i = 0; i < numChares; i++) {
      msgs[i] = new (n/numChares) fftMsg;
      msgs[i]->source = thisIndex;
    }

#ifdef MODE_CPU
    // Reduction to the mainchare to signal that initialization is complete
    contribute(CkCallback(CkReductionTarget(Main, FFTReady), mainProxy));
#endif
  }

  void initComplete() {
#ifdef MODE_CUDA
    contribute(CkCallback(CkReductionTarget(Main, FFTReady), mainProxy));
#endif
  }

  void prepTranspose() {
#ifdef MODE_CPU
    thisProxy[thisIndex].prepTransposeDone();
#elif defined MODE_CUDA
    hapiCheck(cudaEventRecord(compute_event, compute_stream));
    hapiCheck(cudaStreamWaitEvent(comm_stream, compute_event, 0));

    if (iteration == 0) {
      hapiCheck(cudaMemcpyAsync(h_in, d_in, sizeof(complex_t) * n,
            cudaMemcpyDeviceToHost, comm_stream));
    } else {
      hapiCheck(cudaMemcpyAsync(h_out, d_out, sizeof(complex_t) * n,
            cudaMemcpyDeviceToHost, comm_stream));
    }
    CkCallback* cb = new CkCallback(CkIndex_fft::prepTransposeDone(), thisProxy[thisIndex]);
    hapiAddCallback(comm_stream, cb);
#endif
  }

  void sendTranspose() {
#ifdef MODE_CPU
    complex_t* src_buf = (iteration == 0) ? in : out;
#elif defined MODE_CUDA
    complex_t* src_buf = (iteration == 0) ? h_in : h_out;
#endif

    // All-to-all transpose by constructing and sending
    // point-to-point messages to each chare in the array.
    for (int i = thisIndex; i < thisIndex+numChares; i++) {
      // Stagger communication order to avoid hotspots and the
      // associated contention.
      int k = i % numChares;
      for (int j = 0, l = 0; j < N/numChares; j++) {
#ifdef MODE_CPU
        memcpy(msgs[k]->data[(l++)*N/numChares], src_buf[k*N/numChares+j*N],
            sizeof(complex_t)*N/numChares);
#elif defined MODE_CUDA
        memcpy(&msgs[k]->data[(l++)*N/numChares], &src_buf[k*N/numChares+j*N],
            sizeof(complex_t)*N/numChares);
#endif
      }

      // Tag each message with the iteration in which it was
      // generated, to prevent mis-matched messages from chares that
      // got all of their input quickly and moved to the next step.
      CkSetRefNum(msgs[k], iteration);
      thisProxy[k].getTranspose(msgs[k]);
      // Runtime system takes ownership of messages once they're sent
      msgs[k] = NULL;
    }
  }

  void applyTranspose(fftMsg *m) {
    int k = m->source;
#ifdef MODE_CPU
    for (int j = 0, l = 0; j < N/numChares; j++) {
      for (int i = 0; i < N/numChares; i++) {
        out[k*N/numChares+(i*N+j)][0] = m->data[l][0];
        out[k*N/numChares+(i*N+j)][1] = m->data[l++][1];
      }
    }
#elif defined MODE_CUDA
    // Tranpose on CPU
    for (int j = 0, l = 0; j < N/numChares; j++) {
      for (int i = 0; i < N/numChares; i++) {
        h_out[k*N/numChares+(i*N+j)].x = m->data[l].x;
        h_out[k*N/numChares+(i*N+j)].y = m->data[l++].y;
      }
    }

    /*
    // TODO: Transpose on GPU
    // Memcpy received data to pinned host buffer
    memcpy(&h_out[k*n/numChares], m->data[0], sizeof(complex_t) * n/numChares);

    // Transfer to device memory
    hapiCheck(cudaMemcpyAsync(&d_out[k*n/numChares], &h_out[k*n/numChares],
          sizeof(complex_t) * n/numChares, cudaMemcpyHostToDevice, comm_stream));

    // FIXME: Need to be on compute stream?
    invokeTranspose();
    */
#endif

    // Save just-received messages to reuse for later sends, to
    // avoid reallocation
    delete msgs[k];
    msgs[k] = m;
    msgs[k]->source = thisIndex;
  }

  void transferTransposed() {
#ifdef MODE_CUDA
    // Transfer to device memory
    hapiCheck(cudaMemcpyAsync(d_out, h_out, sizeof(complex_t) * n,
          cudaMemcpyHostToDevice, comm_stream));
#endif
  }

  void fftExecute() {
#ifdef MODE_CPU
    fftw_execute(p1);
#elif defined MODE_CUDA
    hapiCheck(cudaEventRecord(comm_event, comm_stream));
    hapiCheck(cudaStreamWaitEvent(compute_stream, comm_event, 0));

    if (cufftExecC2C(p1, d_out, d_out, CUFFT_FORWARD) != CUFFT_SUCCESS) {
      CkAbort("CUFFT Error: Unable to execute plan");
    }
#endif
  }

  void twiddle(double sign) {
    int k = thisIndex;
#ifdef MODE_CPU
    double a, c, s, re, im;
    for (int i = 0; i < N/numChares; i++) {
      for (int j = 0; j < N; j++) {
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
#elif defined MODE_CUDA
    invokeTwiddle(d_out, N, numChares, k, sign, compute_stream);
#endif
  }

  void initValidation() {
    validating = true;

#ifdef MODE_CPU
    memcpy(in, out, sizeof(complex_t) * n);

    fftw_destroy_plan(p1);
    int length[] = {N};
    p1 = fftw_plan_many_dft(1, length, N/numChares, out, length, 1, N,
                            out, length, 1, N, FFTW_BACKWARD, FFTW_ESTIMATE);
#elif defined MODE_CUDA
    // TODO
#endif

    contribute(CkCallback(CkReductionTarget(Main, FFTReady), mainProxy));
  }

  void calcResidual() {
#ifdef MODE_CPU
    double infNorm = 0.0;

    srand48(thisIndex);
    for(int i = 0; i < n; i++) {
      out[i][0] = out[i][0]/(N*N) - drand48();
      out[i][1] = out[i][1]/(N*N) - drand48();

      double mag = sqrt(pow(out[i][0], 2) + pow(out[i][1], 2));
      if(mag > infNorm) infNorm = mag;
    }

    double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));
#elif defined MODE_CUDA
    // TODO
    double r = 0;
#endif

    CkCallback cb(CkReductionTarget(Main, printResidual), mainProxy);
    contribute(sizeof(double), &r, CkReduction::max_double, cb);
  }

  fft(CkMigrateMessage* m) {}
  ~fft() {}
};

#include "fft1d.def.h"
