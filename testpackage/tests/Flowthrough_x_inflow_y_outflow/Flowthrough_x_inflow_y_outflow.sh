#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="bulk.0000001.vlsv"

# Serial execution
SHA_ref="fe41604c502f12f2406af1984cc103bfec36b564f83ea26ea40a3b789fb21095"

source ../../run_and_check.sh $PWD