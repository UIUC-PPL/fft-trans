//#include <stdio.h>
//#include <iostream>

void readCommFile(fftw_complex *data, char *filename)
{
  FILE *pFile;
  if(!(pFile = fopen (filename,"r"))){
    //printf("Warning: File not found or open failure on rank %d\n",rank);
    printf("File open failed\n");
    //MPI_Abort(MPI_COMM_WORLD,1);
    return;
  }

  int l = 0;
  while(fscanf (pFile, "%lf %lf", &data[l][0],&data[l][1]) != EOF) {l++;}

  fclose(pFile);
}

void writeCommFile(int n, fftw_complex *data, char *filename)
{
  FILE *pFile;
  if(!(pFile = fopen (filename,"w"))){
    //printf("Warning: File not found or open failure on rank %d\n",rank);
    printf("File open for write failed\n");
    //MPI_Abort(MPI_COMM_WORLD,1);
    return;
  }

  for(int l=0; l<n; l++)
    fprintf(pFile, "%.24f %.24f\n", data[l][0],data[l][1]);

  fclose(pFile);
}

