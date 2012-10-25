#!/bin/bash

PROCS=$1
N=$2

./charmrun +p${PROCS} ++local main ${N}
