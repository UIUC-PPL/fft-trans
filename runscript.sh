#!/bin/bash

NODES=$1
PROCS=$1

qsub -t 10 -n ${NODES} -A CharmRTS --mode vn ./fft1d $2 $3

#qsub -A PARTS -n 256 -t 15 --mode vn --env BG_MAPPING=ZYXT ./fft1d