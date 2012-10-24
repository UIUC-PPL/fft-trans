#!/bin/bash

NODES=$1
PROCS=$2
NUMBUF=$3

qsub -t 10 -n ${NODES} -A CharmRTS --mode vn --env BG_MAPPING=ZYXT ./fft1d${PROCS} ${NUMBUF}
