#!/bin/sh

MODE=vn
MAPPING=ZYXT
NP=32768
BINARY=fft1d32k

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 416
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 448
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 480
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 512
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 544
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 576
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 608

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 416
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 512
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 608
