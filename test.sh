#!/bin/bash

./charmrun +p$1 ++local fft1d
#mpirun -np $1 ./fft_ref $2
