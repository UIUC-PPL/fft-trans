#!/bin/bash

./charmrun +p1 ++local fft1d $1 $2
mpirun -np $1 ./fft_ref $2
