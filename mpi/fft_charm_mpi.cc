#include <mpi.h> 
#include <cstdlib>
#include <cmath>
#include <limits>
#include <fftw3.h>
#include "fft1d_mpi_wrapper.h"
#include "mpi-interoperate.h"

using namespace std;

int main(int argc, char *argv[]) {
  int rank, size; 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

  uint64_t N = atoi(argv[1]);
  uint64_t local_size=N*N/size;
  int validate = 1; 
  //int validate = atoi(argv[2]);

  CharmLibInit(MPI_COMM_WORLD, argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);

  fftw_complex *in, *out;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * local_size);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * local_size);

  srand48(rank);
  for(uint64_t i = 0; i < local_size; i++) {
    in[i][0] = drand48();
    in[i][1] = drand48();
  }

  //call fft
  fft1d(in, out, N);
  MPI_Barrier(MPI_COMM_WORLD);

  if(validate) {
    double infNorm = 0.0;
    srand48(rank);
    for(uint64_t i = 0; i < local_size; i++) {
      out[i][0] = out[i][0]/(N*N) - drand48();
      out[i][1] = out[i][1]/(N*N) - drand48();

      double mag = sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]);
      if(mag > infNorm) infNorm = mag;
    }

    double my_r = infNorm / (std::numeric_limits<double>::epsilon() * log((double)N * N));
    double r;

    MPI_Reduce(&my_r, &r, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0) {
      if(r < 16)
        printf("r = %g, PASS!\n",r);
      else
        printf("r = %g, FAIL\n",r);
    }
  }

  fftw_free(in);
  fftw_free(out);

  CharmLibExit();
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();      
  return 0; 
} 
