#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="bulk.0000001.vlsv"

# Serial execution
SHA_ref="f2923928ca9a41ee6e97b725818684826e3d2bd94e780bfac3917e040ab15883"

source ../../run_and_check.sh $PWD