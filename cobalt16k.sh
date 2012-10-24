#!/bin/sh

MODE=vn
MAPPING=ZYXT
NP=16384
BINARY=fft1d16k

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 208
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 224
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 240
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 256
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 272
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 288
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 304

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 208
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 256
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 304
