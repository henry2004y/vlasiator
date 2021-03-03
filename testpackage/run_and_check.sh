#!/bin/bash
# To be sourced by each test script.
# ./runtest.sh test_folder_directory

echo "Start testing."
cd "$1"
configfile=$(ls | grep *.cfg)

if [ ! -d ../reference ]
then
   scp -o 'ProxyCommand ssh hongyang@login.physics.helsinki.fi -W %h:%p' \
   hongyang@turso.cs.helsinki.fi:proj/reference.tar.gz ../
   tar -xzf ../reference.tar.gz -C ../
fi

export OMP_NUM_THREADS=$nthreads
mpiexec -n $nprocs ../../../vlasiator --run_config $configfile

isIdentical=$(julia -e ' using Vlasiator
   filenames, tol = ARGS[1:2], parse(Float64, ARGS[3])
   isIdentical = compare(filenames[1], filenames[2], tol)
   println(isIdentical)
   ' ../reference/$(basename $1)/$filecheck $filecheck $tol)

if [ "$isIdentical" = "true" ]
then
   rm *vlsv *txt
else
   exit 1   
fi