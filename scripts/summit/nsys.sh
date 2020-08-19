#!/bin/bash
#BSUB -W 10
#BSUB -P csc357
#BSUB -nnodes 1
#BSUB -J fft1d-nsys

# These need to be changed between submissions
file=fft1d
n_nodes=1
n_procs=$((n_nodes * 6))
N=12288

# Function to display commands
exe() { echo "\$ $@" ; "$@" ; }

cd $HOME/work/fft-trans

ppn=1
pemap="L0,4,8,84,88,92"

echo "# 1D FFT Nsight Profiling"

for mul in 1 2 4 8 16
do
  n_chares=$((n_procs * mul))
  echo "# Multiplier $mul ($n_chares chares on $n_procs processes)"
  exe jsrun -n$n_procs -a1 -c$ppn -g1 -K3 -r6 nsys profile -f true -o fft1d-N$N-m$mul-n$n_nodes-p%q{OMPI_COMM_WORLD_RANK} ./$file -c $n_chares -n $N +ppn $ppn +pemap $pemap
done
