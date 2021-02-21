#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="bulk.0000001.vlsv"

# Serial execution
SHA_ref="19eb05092c85a64e7f70ae3c652e4055eac026bddb6dc1654ed290dca27d0fbc"

source ../../run_and_check.sh $PWD
