#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="fullf.0000001.vlsv"

# Number of MPI processes
nprocs=1
# Number of OpenMP threads
nthreads=1
# Tolerance that scales with quantity magnitude
tol=1e-2

source ../../run_and_check.sh $PWD
