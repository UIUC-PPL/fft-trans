#!/bin/bash

./charmrun +p$1 ++local fft1d $2
mpirun -np $1 ./fft_ref $2
