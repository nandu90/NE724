#!/bin/bash

PROCS=$1
ARG1=$2
ARG2=$3
rm -rf RUN
echo "Running Job on $PROCS processors"
echo "See RUN folder"
mkdir RUN
cd RUN
mpirun -np $PROCS ../bin/system1PH $ARG1 $ARG2


