#!/bin/bash
# Run the test and compare the result.

cd $(dirname "$0")

# Compared filename
filecheck="bulk.0000001.vlsv"

# Serial execution
SHA_ref="99d75cb65f9e996a868e891a19b8f9ee882222307c349bd539d2834f6a1daf1e"

source ../../run_and_check.sh $PWD
