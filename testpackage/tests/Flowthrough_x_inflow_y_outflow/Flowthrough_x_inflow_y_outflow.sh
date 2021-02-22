#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="bulk.0000001.vlsv"

# Serial execution
SHA_ref="4984d41e6772d452276ad095cb5437d23842f5fe0ca6617f9ba46ddeff186437"

source ../../run_and_check.sh $PWD
