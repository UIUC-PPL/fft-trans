#!/bin/bash

NODES=$1
PROCS=$2

qsub -t 10 -n ${NODES} -A CharmRTS --mode c${PROCS} -q R.l2p ./fft1d $3 $4
