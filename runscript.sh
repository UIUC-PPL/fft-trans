#!/bin/bash

NODES=$1
PROCS=$1

qsub -t 10 -n ${NODES} -A CharmRTS --mode vn ./fft1d $2 $3
