#include <fftw3.h>
#include "NDMeshStreamer.h"

#define TWOPI 6.283185307179586
#define SET_VALUES(a,b,c)  do { (a)[0] = b; (a)[1] = c; } while (0);

#include "fft.decl.h"

struct fft : public MeshStreamerGroupClient<fftw_complex> {
  fft_SDAG_CODE

  int iteration, count, sign;
  uint64_t n, N;
  fftw_plan p1;
  fftw_complex *in, *out, *buf;
  streamer_t streamer;

  fft(uint64_t N, CProxy_fftData data, int sign, CkCallback startCB) : n(N/CkNumPes()), N(N), sign(sign) {
    in = data.ckLocalBranch()->getIn();
    out = data.ckLocalBranch()->getOut();

    p1 = fftw_plan_many_dft(1, (int*)&N, n, out, (int*)&N, 1, N,
                            out, (int*)&N, 1, N, sign, FFTW_ESTIMATE);

    buf = new fftw_complex[n*n];

    // Reduction to signal that initialization is complete
    contribute(startCB);
  }

  void initStreamer(streamer_t streamer_) {
    streamer = streamer_;
    streamer.ckLocalBranch()->init(1, CkCallback(CkIndex_fft::streamerReady(), thisProxy), CkCallback(CkIndex_fft::doneStreaming(), thisProxy), 0, false);
  }

  void sendTranspose(fftw_complex *src_buf) {
    for(int i = 0; i < CkNumPes(); i++) {
      for(int j = 0, l = 0; j < n; j++)
        memcpy(buf[(l++)*n], src_buf[i*n+j*N], sizeof(fftw_complex)*n);

      streamer.ckLocalBranch()->insertData(buf, n*n, i);
    }
    streamer.ckLocalBranch()->done();
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
};

#include "fft.def.h"
