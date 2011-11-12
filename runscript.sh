#!/bin/bash

NODES=$1
PROCS=$2

qsub -t 10 -n ${NODES} -A CharmRTS --mode vn ./fft1d${PROCS}
