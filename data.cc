#include <limits>

struct fftData : public CBase_fftData {
  fftw_complex *in, *out;
  uint64_t n, N;

  fftData(uint64_t N) : n(N/CkNumPes()), N(N) {
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n*N);
    srand48(CkMyPe());
    for(int i = 0; i < n*N; i++) {
      in[i][0] = drand48();
      in[i][1] = drand48();
    }
  }

  fftw_complex* getIn() { return in; }
  fftw_complex* getOut() { return out; }

  void swap(CkCallback cb) { memcpy(in, out, sizeof(fftw_complex) * n*N); contribute(cb); }

  void calcResidual(CkCallback cb) {
    double infNorm = 0.0;

    srand48(CkMyPe());
    for(int i = 0; i < n*N; i++) {
      out[i][0] = out[i][0]/(N*N) - drand48();
      out[i][1] = out[i][1]/(N*N) - drand48();

      double mag = sqrt(pow(out[i][0], 2) + pow(out[i][1], 2));
      if(mag > infNorm) infNorm = mag;
    }

    double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));

    contribute(sizeof(double), &r, CkReduction::max_double, cb);
  }
};
