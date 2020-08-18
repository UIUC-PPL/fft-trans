#include "fft1d.h"

void readCommFile(complex_t *data, char *filename) {
  FILE *pFile;
  if (!(pFile = fopen(filename,"r"))) {
    printf("File open failed\n");
    return;
  }

  int l = 0;
#ifdef MODE_CPU
  while (fscanf(pFile, "%lf %lf", &data[l][0], &data[l][1]) != EOF) {l++;}
#elif defined MODE_CUDA
  while (fscanf(pFile, "%f %f", &data[l].x, &data[l].y) != EOF) {l++;}
#endif

  fclose(pFile);
}

void writeCommFile(int n, complex_t *data, char *filename) {
  FILE *pFile;
  if (!(pFile = fopen(filename, "w"))) {
    printf("File open for write failed\n");
    return;
  }

  for (int l = 0; l < n; l++) {
#ifdef MODE_CPU
    fprintf(pFile, "%.24f %.24f\n", data[l][0], data[l][1]);
#elif defined MODE_CUDA
    fprintf(pFile, "%.24f %.24f\n", data[l].x, data[l].y);
#endif
  }

  fclose(pFile);
}
