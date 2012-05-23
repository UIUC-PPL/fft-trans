#!/bin/bash

NODES=$1
PROCS=$1

qsub -t 10 -n ${NODES} -A CharmRTS --mode c8 -q R.l2p ./fft1d $2 $3
