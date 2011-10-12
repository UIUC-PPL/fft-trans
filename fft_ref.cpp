#include <stdio.h>
#include <iostream>

#include <mpi.h> 
#include <cstdlib>
#include "fftw3-mpi.h"

using namespace std;

void readCommFile(char *filename);

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

  ptrdiff_t alloc_local =  fftw_mpi_local_size_1d(N*N, MPI_COMM_WORLD,
      FFTW_FORWARD, FFTW_ESTIMATE, &local_ni, &local_i_start,
      &local_no, &local_o_start);

  data = fftw_alloc_complex(alloc_local);

  plan = fftw_mpi_plan_dft_1d(N*N, data, data, MPI_COMM_WORLD,FFTW_FORWARD, FFTW_ESTIMATE);

  for (int i=0; i<N*N/size; i++) {
    data[i][0] = rank*N*N/size+i;
    data[i][1] = 0.0;
    printf("init: [%d].%d = %f\n",rank,i,data[i][0]);
  }

  /*
  char filename[80];
  int x2,y2,z2,send_size,recv_size;
  FILE * pFile;

  int *send_sizes = new int[size];

  sprintf(filename,"%d-%d.dump%d",size,N,rank);
  readCommFile(filename);
  */

  fftw_execute(plan);

  for (int i=0; i<N*N/size; i++)
    printf("[%d] in[%d] = %f + %fi\n",rank, i, data[i][0], data[i][1]);

  fftw_destroy_plan(plan);

  fftw_mpi_cleanup();
  MPI_Finalize();      
  return 0; 
} 

void readCommFile(char *filename)
{
  FILE *pFile;
  if(!(pFile = fopen (filename,"r"))){
    //printf("Warning: File not found or open failure on rank %d\n",rank);
    //MPI_Abort(MPI_COMM_WORLD,1);
    return;
  }

  int x2, z2, y2, comm_volume;
  unsigned int rank;

  while(fscanf (pFile, "%d %d %d %d", &x2,&y2,&z2,&comm_volume) != EOF)
  {
  }

  fclose(pFile);
}
