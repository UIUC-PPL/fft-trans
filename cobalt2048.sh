#!/bin/sh

MODE=vn
MAPPING=ZYXT
NP=2048
BINARY=fft1d2048

cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 29
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 30
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 31
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 32
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 33
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 34
cobalt-mpirun -mode ${MODE} -env BG_MAPPING=${MAPPING} -np ${NP} ${BINARY} 35
