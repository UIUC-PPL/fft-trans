#include <limits>

void fftData::calcResidual(CkCallback cb) {
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
