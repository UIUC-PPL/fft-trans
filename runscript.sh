#!/bin/bash

NODES=$1
PROCS=$2

qsub -t 10 -n ${NODES} -A CharmRTS --mode c${PROCS} ./fft1d $3 $4
