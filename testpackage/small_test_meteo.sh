#!/bin/bash -l
#PBS -l mppwidth=24
#PBS -l mppnppn=12
#PBS -l walltime=00:15:00
#PBS -V  
#PBS -N test_small


#threads
t=12  

#command for running stuff
run_command="aprun -n 2 -N 1 -d $t"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#folder for reference data 
reference_dir="/stornext/field/users/hoilijo/vlasiator/reference_data/revision1956"


#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh