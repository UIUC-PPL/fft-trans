#include <limits>

void fft::initValidation() {
  memcpy(in, out, sizeof(fftw_complex) * n*N);

  validating = true;
  fftw_destroy_plan(p1);
  p1 = fftw_plan_many_dft(1, (int*)&N, n, out, (int*)&N, 1, N,
                          out, (int*)&N, 1, N, FFTW_BACKWARD, FFTW_ESTIMATE);

  doFFT(CkCallback(CkCallback::ignore));
}

void fft::calcResidual() {
  double infNorm = 0.0;

  srand48(CkMyPe());
  for(int i = 0; i < n*N; i++) {
    out[i][0] = out[i][0]/(N*N) - drand48();
    out[i][1] = out[i][1]/(N*N) - drand48();

    double mag = sqrt(pow(out[i][0], 2) + pow(out[i][1], 2));
    if(mag > infNorm) infNorm = mag;
  }

  double r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));

  CkCallback cb(CkReductionTarget(fft, printResidual), thisProxy[0]);
  contribute(sizeof(double), &r, CkReduction::max_double, cb);
}

void fft::printResidual(double r) {
  CkPrintf("residual = %g\n", r);
  CkExit();
}
