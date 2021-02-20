#!/bin/bash

export OMP_NUM_THREADS=1
export MPICH_MAX_THREAD_SAFETY=funneled

#umask 007

run="mpiexec -n 1"
run_tools="mpiexec -n 1"

# If 1, the reference vlsv files are generated;
# otherwise, compare with referenced files
create_verification_files=0

# folder for all reference data 
reference_dir="ref"
# Compare agains which revision in the first 6 chars of SHA
reference_revision="a2b6b5"

# Define test ------------------------------------------------------------------

## Define test and runs

exec=$( readlink -f ../vlasiator )
diff=$( readlink -f ../vlsvdiff_DP )

# where the tests are run
run_dir="run"

# where the directories for different tests, including cfg and other needed data
# files are located 
test_dir="tests"

# Flowthrough tests
test_name="Flowthrough_trans_periodic"
comparison_vlsv="bulk.0000001.vlsv"
comparison_phiprof="phiprof_0.txt"
variable_names="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components="0 0 1 2 0 1 2 0 1 2"

# Run tests --------------------------------------------------------------------
## add absolute paths to folder names, filenames
reference_dir=$( readlink -f $reference_dir )
run_dir=$( readlink -f $run_dir )_$( date +%Y.%m.%d_%H.%M.%S)
test_dir=$( readlink -f $test_dir)


flags=$( $run $exec  --version | grep CXXFLAGS)
solveropts=$(echo $flags|sed 's/[-+]//g' | gawk '{for(i = 1;i<=NF;i++) { if( $i=="DDP" || $i=="DFP" || index($i,"PF")|| index($i,"DVEC") || index($i,"SEMILAG") ) printf "__%s", $(i) }}')
revision=$( $run $exec --version | gawk '{if(flag==1) {print $1;flag=0}if ($3=="log") flag=1;}' | head -c 6 )

if [ $create_verification_files == 1 ]
then
    # If we create the references, then lets simply run in the reference dir and
    # turn off tests below. Revision is 
    # automatically obtained from the --version output
    reference_revision=${revision}${solveropts}
    echo "Computing reference results into ${reference_dir}/${reference_revision}"
fi

mkdir -p $run_dir 

echo "running ${test_name}"
# directory for test results
vlsv_dir=${run_dir}/${test_name[$i]}
cfg_dir=${test_dir}/${test_name[$i]}
    
# Check if folder for new run exists, if not create them, otherwise delete old results
if [ ! -d ${vlsv_dir} ]; then
    mkdir -p ${vlsv_dir}
else
    rm -f ${vlsv_dir}/*
fi
    
# change to run directory of the test case
cd ${vlsv_dir}
cp ${cfg_dir}/* .

# Run prerequisite script, if it exists
test -e test_prelude.sh && ./test_prelude.sh

# Run the actual simulation
$run $exec --version  > VERSION.txt
$run $exec --run_config=${test_name[$i]}.cfg

# Run postprocessing script, if it exists
test -e test_postproc.sh && ./test_postproc.sh

# Copy new reference data to correct folder
if [ $create_verification_files == 1 ]
then
    result_dir=${reference_dir}/${reference_revision}/${test_name[$i]}
    if [ -e  $result_dir ]
    then
        echo "remove old results"
        rm -rf $result_dir
    fi

    mkdir -p $result_dir
    cp * $result_dir      
else
    # Compare test case with reference solutions
    echo "--------------------------------------------------------------------------------------------" 
    echo "${test_name}  -  Verifying ${revision}_$solveropts against $reference_revision"    
    echo "--------------------------------------------------------------------------------------------" 
    result_dir=${reference_dir}/${reference_revision}/${test_name}

    #print header
    echo "------------------------------------------------------------"
    echo " ref-time     |   new-time       |  speedup                |"
    echo "------------------------------------------------------------"
	 if [ -e  ${result_dir}/${comparison_phiprof} ] 
	 then
        refPerf=$(grep "Propagate   " ${result_dir}/${comparison_phiprof} | gawk  '(NR==1){print $11}')
	 else
	     refPerf="NA"
	 fi
	 if [ -e ${vlsv_dir}/${comparison_phiprof} ] 
	 then
        newPerf=$(grep "Propagate   " ${vlsv_dir}/${comparison_phiprof} | gawk  '(NR==1){print $11}')
	 else
	     newPerf="NA"
	 fi
	 # Print speedup if both refPerf and newPerf are numerical values
    speedup=$( echo $refPerf $newPerf |gawk '{if($2 == $2 + 0 && $1 == $1 + 0 ) print $1/$2; else print "NA"}')
    echo  "$refPerf        $newPerf         $speedup"
    echo "------------------------------------------------------------"
    echo "  variable     |     absolute diff     |     relative diff | "
    echo "------------------------------------------------------------"

	 variables=(${variable_names// / })
    #variables=($variable_names[$i])
    #read -a variables <<< variable_names[$i]
	 indices=(${variable_components// / })

    for i in ${!variables[*]}
    do
        if [ "${variables[$i]}" == "fg_e" ] || [ "${variables[$i]}" == "fg_b" ]
        then
            relativeValue=$($diff --meshname=fsgrid  ${result_dir}/${comparison_vlsv[$i]} ${vlsv_dir}/${comparison_vlsv[$i]} ${variables[$i]} ${indices[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
            absoluteValue=$($diff --meshname=fsgrid  ${result_dir}/${comparison_vlsv[$i]} ${vlsv_dir}/${comparison_vlsv[$i]} ${variables[$i]} ${indices[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
            # Print the results
            echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
            
        elif [ ! "${variables[$i]}" == "proton" ]
        then
            relativeValue=$($diff ${result_dir}/${comparison_vlsv[$i]} ${result_dir}/${comparison_vlsv[$i]} \
                ${variables[$i]} ${indices[$i]} | \
                grep "The absolute 0-distance between both datasets" | \
                gawk '{print $8}'  )
            absoluteValue=$($diff ${result_dir}/${comparison_vlsv[$i]} ${vlsv_dir}/${comparison_vlsv[$i]} \
                ${variables[$i]} ${indices[$i]} | \
                grep "The absolute 0-distance between both datasets" | \
                gawk '{print $8}'  )
            # Print the results
            echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
        elif [ "${variables[$i]}" == "proton" ]
        then
            echo "--------------------------------------------------------------------------------------------" 
            echo "   Distribution function diff                                                               "
            echo "--------------------------------------------------------------------------------------------" 
            $diff ${result_dir}/${comparison_vlsv[$i]} ${vlsv_dir}/${comparison_vlsv[$i]} proton 0
        fi
    done # loop over variables

    echo "--------------------------------------------------------------------------------------------" 
fi
