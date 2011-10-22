#!/bin/bash

NODES=1
PROCS=1

qsub -t 10 -n ${NODES} -A CharmRTS --mode smp ./fft_bench $1
