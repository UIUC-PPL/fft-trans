#!/bin/sh

MODE=vn
MAPPING=ZYXT
NP=4096
BINARY=fft1d4096

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 58
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 60
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 62
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 64
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 66
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 68
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 70

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 58
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 64
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY}-mpi 70
