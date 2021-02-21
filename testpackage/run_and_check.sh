#!/bin/bash
# To be sourced by each test script.
# ./runtest.sh test_folder_directory

echo "Start testing."
cd "$1"
configfile=$(ls | grep *.cfg)

mpiexec -n 1 ../../../vlasiator --run_config $configfile

SHA_run=$(cat "$filecheck" | sha256sum | head -c 64)

if [ ! "$SHA_run" = "$SHA_ref" ]
then
   echo "There are differences in the output!"; exit 1
else
   rm *vlsv *txt
fi
