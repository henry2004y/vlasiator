#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="bulk.0000001.vlsv"

# Serial execution
SHA_ref="97033a661653522a46c9d97f288f54a0869c256eb5b5d762c97b69f8d0d01a12"

source ../../run_and_check.sh $PWD
