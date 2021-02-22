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

mpiexec -n 1 ../../../vlasiator --run_config $configfile

isIdentical=$(julia -e '
   using Vlasiator
   filenames = ARGS
   isIdentical = compare(filenames[1], filenames[2])
   println(isIdentical)
   ' ../reference/$1/$filecheck $filecheck)

if [ "$isIdentical" = "false" ]
then
   echo "There are differences in the output!"; exit 1
else
   rm *vlsv *txt
fi

#SHA_run=$(cat "$filecheck" | sha256sum | head -c 64)

#if [ ! "$SHA_run" = "$SHA_ref" ]
#then
#   echo "There are differences in the output!"; exit 1
#else
#   rm *vlsv *txt
#fi
