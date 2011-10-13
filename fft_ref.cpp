//#include <stdio.h>
//#include <iostream>

#include <mpi.h> 
#include <cstdlib>
#include "fftw3-mpi.h"
#include <cmath>
#include <limits>

#include "verify.h"

using namespace std;

//void readCommFile(fftw_complex *data, char *filename);

int main(int argc, char *argv[]) 
{ 
  int rank, size; 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

  fftw_plan plan;
  fftw_complex *data;

  fftw_mpi_init();

  if(rank==0){
    if(argc < 2){
      printf("Usage: ./binary <path to files>\n");
      MPI_Abort(MPI_COMM_WORLD,-1);
    }
  }

  int N = atoi(argv[1]);

  ptrdiff_t local_ni=N*N/size, local_i_start = N*N/size*rank;
  ptrdiff_t local_no=local_ni, local_o_start = local_i_start;

  int b_or_f = FFTW_BACKWARD;

  ptrdiff_t alloc_local =  fftw_mpi_local_size_1d(N*N, MPI_COMM_WORLD,
      b_or_f, FFTW_ESTIMATE, &local_ni, &local_i_start,
      &local_no, &local_o_start);

  data = fftw_alloc_complex(alloc_local);

  plan = fftw_mpi_plan_dft_1d(N*N, data, data, MPI_COMM_WORLD,b_or_f, FFTW_ESTIMATE);

  char filename[80];
  sprintf(filename,"%d-%d.dump%d",size,N,rank);
  readCommFile(data, filename);

  for (int i=0; i<N*N/size; i++) {
    //data[i][0] = rank*N*N/size+i;
    //data[i][1] = 0.0;
    printf("init: [%d].%d = %f + %fi\n",rank,i,data[i][0], data[i][1]);
  }

  fftw_execute(plan);

  for (int i=0; i<N*N/size; i++)
    printf("[%d] in[%d] = %f + %fi\n",rank, i, data[i][0], data[i][1]);
  double infNorm = 0.0;
  srand48(rank);
  for (int i=0; i<N*N/size; i++){
    data[i][0] = data[i][0]/(N*N) - drand48(); 
    data[i][1] = data[i][1]/(N*N) - drand48();

    if(fabs(data[i][0]) > infNorm)
      infNorm = fabs(data[i][0]);
    if(fabs(data[i][1]) > infNorm)
      infNorm = fabs(data[i][1]);

    //printf("[%d] in[%d] = %g + %gi\n",rank, i, data[i][0], data[i][1]);
  }

  double r = infNorm/(std::numeric_limits<double>::epsilon()*log(N*N));

  if(r < 16)
    printf("r = %g, PASS!\n",r);
  else
    printf("r = %g, FAIL\n",r);

  fftw_destroy_plan(plan);

  fftw_mpi_cleanup();
  MPI_Finalize();      
  return 0; 
} 

/*
void readCommFile(fftw_complex *data, char *filename)
{
  FILE *pFile;
  if(!(pFile = fopen (filename,"r"))){
    //printf("Warning: File not found or open failure on rank %d\n",rank);
    printf("File open failed\n");
    MPI_Abort(MPI_COMM_WORLD,1);
    return;
  }

  int l = 0;
  while(fscanf (pFile, "%lf %lf", &data[l][0],&data[l][1]) != EOF) {l++;}

  fclose(pFile);
}
*/
