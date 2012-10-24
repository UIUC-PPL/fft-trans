#!/bin/sh

MODE=vn
MAPPING=ZYXT
NP=1024
BINARY=fft1d1024

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 13
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 14
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 15
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 16
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 17
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 18
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 19
