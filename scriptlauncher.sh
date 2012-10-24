#!/bin/bash

NODES=$1
PROCS=$2

qsub -t 60 -n ${NODES} -A PEACEndStation --mode script cobalt${PROCS}.sh
