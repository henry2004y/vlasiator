#!/usr/bin/env python
# ------------------------------------------------------------------------------
# configure.py: Vlasiator configuration script in python.
#
# It uses the command line options and default settings to create customized
# versions of Makefile from the template.
# Original version by CJW. Modified by Hongyang Zhou.
# ------------------------------------------------------------------------------

# Modules
import argparse
import glob
import re
import subprocess
import os
import pkg_resources
import sys
import warnings
import fileinput

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

# --- Step 1. Prepare parser, add each of the arguments ------------------
vlasiator_description = (
    "Prepare custom Makefile for compiling Vlasiator"
)
vlasiator_epilog = (
    "Full documentation of options available at "
    "https://github.com/fmihpc/vlasiator/wiki/Installing-Vlasiator"
)
parser = argparse.ArgumentParser(description=vlasiator_description, epilog=vlasiator_epilog)

parser.add_argument('-install',
                    action='store_true',
                    default=False,
                    help='install all the dependencies into lib')

parser.add_argument('--nx',
                    type=int,
                    default=4,
                    help='set number of cells in the block x dimension')

parser.add_argument('--ny',
                    type=int,
                    default=4,
                    help='set number of cells in the block y dimension')

parser.add_argument('--nz',
                    type=int,
                    default=4,
                    help='set number of cells in the block z dimension')

parser.add_argument('--vorder',
                    type=int,
                    default=5,
                    choices=[2,3,5],
                    help='set order of velocity space acceleration')

parser.add_argument('--sorder',
                    type=int,
                    default=3,
                    choices=[2,3,5],
                    help='set order of spatial translation')

parser.add_argument('--fieldsolver',
                    default='londrillo_delzanna',
                    help='select field solver')

parser.add_argument('--fieldorder',
                    type=int,
                    default=2,
                    choices=[1,2],
                    help='select field solver order')

parser.add_argument('-amr',
                    action='store_true',
                    default=False,
                    help='enable AMR in velocity space')

parser.add_argument('-noamr', dest='amr', action='store_false')

parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable debug flags; override other compiler options')

parser.add_argument('-nodebug', dest='debug', action='store_false')

parser.add_argument('-debugsolver',
                    action='store_true',
                    default=False,
                    help='enable debug preprocessor for field solver')

parser.add_argument('-nodebugsovler', dest='debugsolver', action='store_false')

parser.add_argument('-debugionosphere',
                    action='store_true',
                    default=False,
                    help='enable debug preprocessor for ionosphere module')

parser.add_argument('-nodebugionosphere', dest='debugionosphere', action='store_false')

parser.add_argument('-debugfloat',
                    action='store_true',
                    default=False,
                    help='enable preprocessor for floating point exceptions')

parser.add_argument('-nodebugfloat', dest='debugfloat', action='store_false')

parser.add_argument('-float',
                    action='store_true',
                    default=False,
                    help='enable single precision')

parser.add_argument('-double', dest='float', action='store_false')

parser.add_argument('-distdouble',
                    action='store_true',
                    default=False,
                    help='enable double precision for distribution function')

parser.add_argument('-distsingle', dest='distdouble', action='store_false')

parser.add_argument('-profile', 
                    dest='profile', 
                    action='store_true',
                    default=True,
                    help='enable profiler')

parser.add_argument('-noprofile', dest='profile', action='store_false')

parser.add_argument('-mpi',
                    action='store_true',
                    default=True,
                    help='enable parallelization with MPI')

parser.add_argument('-nompi', dest='mpi', action='store_false')

parser.add_argument('-omp',
                    action='store_true',
                    default=False,
                    help='enable parallelization with OpenMP')

parser.add_argument('-noomp', dest='omp', action='store_false')

parser.add_argument('--machine',
                    default='default',
                    help='choose machine name for specific setup')

parser.add_argument('-papi',
                    action='store_true',
                    default=False,
                    help='enable Papi memory profiler')

parser.add_argument('-jemalloc',
                    action='store_true',
                    default=False,
                    help='enable JEMALLOC allocator')

parser.add_argument('-silo',
                    action='store_true',
                    default=False,
                    help='enable silo format converter')

cxx_choices = [
    'g++',
    'icpc',
    'icpc-debug',
    'cray',
    'clang++',
    'clang++-apple',
]


parser.add_argument(
    '--cxx',
    default='g++',
    choices=cxx_choices,
    help='select C++ compiler and default set of flags')

parser.add_argument(
    '--cflag',
    default=None,
    help='additional string of flags to append to compiler/linker calls')

# Parse command-line inputs
args = vars(parser.parse_args())

# --------- Step 2. Test for incompatible arguments ----------------------------

if args['install'] and args['machine'] != 'default':
    raise SystemExit('### CONFIGURE ERROR: does not support fresh installation with a preset machine makefile')

if not os.path.isfile(os.path.join("MAKE/Makefile."+args['machine'])):
    raise SystemExit('Target machine does not exist!')

# --------- Step 3. Set Makefile options based on above argument ---------------

def display_settings():
    print('Vlasiator has now been configured with the following options:')
    print('  Machine:                    ' + (args['machine'] if args['machine'] else 'new'))
    print('  Floating-point precision:   ' + ('single' if args['float'] else 'double'))
    print('  Distribution precision:     ' + ('double' if args['distdouble'] else 'single'))
    print('  Block size:                 ' + str(args['nx']) + ' ' \
                                           + str(args['ny']) + ' ' \
                                           + str(args['nz']))
    print('  MPI parallelism:            ' + ('ON' if args['mpi'] else 'OFF'))
    print('  OpenMP parallelism:         ' + ('ON' if args['omp'] else 'OFF'))
    print('  Order of field solver:      ' + str(args['fieldorder']))
    print('  Order of semilag velocity:  ' + str(args['vorder']))
    print('  Order of semilag spatial:   ' + str(args['sorder']))
    print('  AMR:                        ' + ('ON' if args['amr'] else 'OFF'))
    print('  Profiler:                   ' + ('ON' if args['profile'] else 'OFF'))
    print('  Memory tracker:             ' + ('ON' if args['papi'] else 'OFF'))
    print('  Debug flags:                ' + ('ON' if args['debug'] else 'OFF'))
    print('  Compiler:                   ' + args['cxx'])


if len(sys.argv) == 1:
    # Check existing Makefile if no argument is passed
    if not os.path.isfile("Makefile"):
        print("Vlasiator is not installed.")
    else:
        with open('Makefile') as f:
            fstr = f.read()[110:400]
            args['machine'] = fstr.split("Makefile.",1)[1].split("\n",1)[0]
            if 'mpicxx' in fstr:
                args['cxx'] = 'mpicxx'
            
            args['omp'] = True if 'openmp' in fstr or 'omp' in fstr else False 
            precision = 'DP' if '-DDP ' in fstr else 'SP'
            args['distdouble'] = True if '-DSPF' in fstr else '-DDPF'
            args['profile'] = True if '-DPROFILE' in fstr else False
            if 'ACC_SEMILAG_PQM' in fstr:
                args['vorder'] = 5
            elif 'ACC_SEMILAG_PPM' in fstr:
                args['vorder'] = 3
            elif 'ACC_SEMILAG_PLM' in fstr:
                args['vorder'] = 2
            
            if 'TRANS_SEMILAG_PQM' in fstr:
                args['sorder'] = 5
            elif 'TRANS_SEMILAG_PPM' in fstr:
                args['sorder'] = 3
            elif 'TRANS_SEMILAG_PLM' in fstr:
                args['sorder'] = 2

            args['fieldorder'] = 1 if '1ST_ORDER_SPACE' in fstr else 2

            args['debug'] = True if 'NDEBUG' in fstr else False

            args['amr'] = True if '-DAMR' in fstr else False
            args['jemalloc'] = True if '-DUSE_JEMALLOC' in fstr else False
            args['papi'] = True if '-DPAPI_MEM' in fstr else False

            display_settings()
            sys.exit()



# Prepare dictionaries of substitutions to be made
makefile_options = {}
makefile_options['PREPROCESSOR_FLAGS'] = ''

if args['cxx'] == 'g++':
    makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['COMPILER_FLAGS'] = (
        '-O3 -fopenmp -funroll-loops '
        '-std=c++17 -W -Wall -Wno-unused' )
if args['cxx'] == 'icpc':
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++17 -ipo -xhost -inline-forceinline -qopenmp-simd '
      '-qopt-prefetch=4 '
      '-qoverride-limits' )
if args['cxx'] == 'icpc-debug':
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -xhost -qopenmp-simd -fp-model precise '
      '-qopt-prefetch=4 '
      '-qopt-report=5 -qopt-report-phase=openmp,vec -g -qoverride-limits'
    )
if args['cxx'] == 'cray':
    makefile_options['COMPILER_COMMAND'] = 'CC'
    makefile_options['COMPILER_FLAGS'] = (
        '-O3 -h std=c++17 -h aggress -h vector3 -hfp3' )
    makefile_options['LINKER_FLAGS'] = '-hwp -hpl=obj/lib'
if args['cxx'] == 'clang++':
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++17'
if args['cxx'] == 'clang++-simd':
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++17 -fopenmp-simd'
if args['cxx'] == 'clang++-apple':
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++17'

if args['float']:
    precision = 'SP'
else:
    precision = 'DP'

makefile_options['PREPROCESSOR_FLAGS'] += ' -D'+precision

for key in ('vlsvdiff','vlsvextract','vlsv2silo'):
    makefile_options[key.upper()] = key+'_'+precision

# Distribution function precision
if args['distdouble']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4D_FALLBACK' # vector backend type
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DSPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4F_FALLBACK' # vector backend type

if args['profile']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DPROFILE'
 
if args['vorder'] == 2:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PLM'       
elif args['vorder'] == 3:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PPM'
elif args['vorder'] == 5:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PQM' 
else:
    raise SystemExit('### CONFIGURE ERROR: unknown semilag solver order for velocity space acceleration')

if args['sorder'] == 2:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PLM'       
elif args['sorder'] == 3:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PPM'
elif args['sorder'] == 5:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PQM' 
else:
    raise SystemExit('### CONFIGURE ERROR: unknown semilag solver order for spatial translation')


if args['fieldorder'] == 1:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DFS_1ST_ORDER_SPACE'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DFS_1ST_ORDER_TIME'

if args['amr']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DAMR'

if args['debug']:
    pass
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DNDEBUG'

if args['debugsolver']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDEBUG_SOLVERS'

if args['debugionosphere']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDEBUG_IONOSPHERE'

# Catch floating point exceptions and stop execution
if args['debugfloat']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DCATCH_FPE'

if args['papi']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DPAPI_MEM'

if args['jemalloc']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DUSE_JEMALLOC -DJEMALLOC_NO_DEMANGLE'

if args['debug']:
    if (args['cxx'] == 'g++' or args['cxx'] == 'icpc' 
        or args['cxx'] == 'icpc-debug' or args['cxx'] == 'clang++'
        or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g'  # -Og
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] = '-O0 -h std=c++11'

if args['mpi']:
    if (args['cxx'] == 'g++' or args['cxx'] == 'icpc' 
        or args['cxx'] == 'icpc-debug' or args['cxx'] == 'clang++'
        or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_COMMAND'] = 'mpicxx'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -h mpi1'
else:
    raise SystemExit('### CONFIGURE ERROR: -mpi is required for compilation!')

if args['omp']:
    if args['cxx'] == 'g++' or args['cxx'] == 'clang++':
        makefile_options['COMPILER_FLAGS'] += ' -fopenmp'
    if args['cxx'] == 'clang++-apple':
        # Apple Clang disables the front end OpenMP driver interface; enable it
        # via the preprocessor. 
        # Must install LLVM's OpenMP runtime library libomp beforehand
        makefile_options['COMPILER_FLAGS'] += ' -Xpreprocessor -fopenmp'
    if args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug':
        makefile_options['COMPILER_FLAGS'] += ' -qopenmp'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -homp'

if args['cflag'] is not None:
    makefile_options['COMPILER_FLAGS'] += ' '+args['cflag']

# Install dependencies if required
if args['install']:
    if not os.path.isdir("lib"):
        subprocess.check_call(["mkdir", "lib"]) 
    if not os.path.isdir("lib/vectorclass"):
        subprocess.check_call(["git", "clone", "https://github.com/vectorclass/version1.git"])
        subprocess.check_call(["git", "clone", "https://github.com/vectorclass/add-on.git"])
        subprocess.check_call(["cp", "add-on/vector3d/vector3d.h", "version1/"])
        subprocess.check_call(["mv", "version1", "lib/vectorclass"])
    if not os.path.isdir("lib/fsgrid"):
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/fsgrid.git"])
        subprocess.check_call(["mv", "fsgrid", "lib"])
    if not os.path.isdir("lib/dccrg"): 
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/dccrg.git"])
        subprocess.check_call(["mv", "dccrg", "lib"])
        os.chdir("lib/dccrg")
        subprocess.check_call(["git", "checkout", "01482cfba8"])
        os.chdir("../..")
    if not os.path.isdir("lib/phiprof"):
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/phiprof.git"])
        subprocess.check_call(["mv", "phiprof", "lib"])
        os.chdir("lib/phiprof/src")
        subprocess.check_call(["make"])
        os.chdir("../../..")
    if not os.path.isdir("lib/vlsv"):
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/vlsv.git"])
        subprocess.check_call(["mv", "vlsv", "lib"])
        os.chdir("lib/vlsv")
        subprocess.check_call(["make"])
        os.chdir("../..")
    if not os.path.isdir("lib/jemalloc"):
        subprocess.check_call(["wget", \
            "https://github.com/jemalloc/jemalloc/releases/download/4.0.4/jemalloc-4.0.4.tar.bz2"])
        subprocess.check_call(["tar", "-xf", "jemalloc-4.0.4.tar.bz2"])
        os.chdir("jemalloc-4.0.4")
        subprocess.check_call(["./configure","--prefix="+os.getcwd()+"/jemalloc", \
        "--with-jemalloc-prefix=je_"])
        subprocess.check_call(["make"])
        subprocess.check_call(["make", "install"])
        subprocess.check_call(["mv", "jemalloc", "../lib"])
        os.chdir("..")
    if not os.path.isdir("lib/Eigen"):
        subprocess.check_call(["wget", \
        "https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.bz2"])
        subprocess.check_call(["tar", "-xf", "eigen-3.2.8.tar.bz2"])
        subprocess.check_call(["cp", "-r", "eigen-3.2.8/Eigen", "lib"])
    if not os.path.isdir("lib/zoltan"):
        subprocess.check_call(["wget", \
        "http://cs.sandia.gov/Zoltan/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz"])
        subprocess.check_call(["tar", "-xf", "zoltan_distrib_v3.83.tar.gz"])
        subprocess.check_call(["mkdir", "lib/zoltan"])
        os.chdir("lib/zoltan")
        subprocess.check_call(["../../Zoltan_v3.83/configure","--prefix="+os.getcwd(), \
        "--enable-mpi", "--with-mpi-compilers", "--with-gnumake", \
        "--with-id-type=ullong"])
        subprocess.check_call(["make", "-j", "4"])
        subprocess.check_call(["make", "install"])
        os.chdir("../..")

    for f in ["add-on", "eigen-3.2.8.tar.bz2","eigen-3.2.8",\
        "zoltan_distrib_v3.83.tar.gz", "Zoltan_v3.83", \
        "jemalloc-4.0.4", "jemalloc-4.0.4.tar.bz2"]:
        subprocess.check_call(["rm", "-rf", f])

    # Boost is skipped as it is too large to install here
    if not os.path.isfile(os.path.join("/usr/lib/x86_64-linux-gnu/libboost_program_options.a")):
        raise SystemExit('Boost not found: try to set the correct path, or install it manually as follows:'
            '\n sudo apt update\n sudo apt install libboost-all-dev')

# --------- Step 4. Create new files, finish up --------------------------------

# Read templates
with open('MAKE/Makefile.in', 'r') as f:
    makefile_template = f.read()
    # Substitute machine name
    makefile_template = makefile_template.replace('@machine@', args['machine'])

    # Make substitutions
    for key, val in makefile_options.items():
        makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

    # Redirect field solver folder
    if args['amr']:
        makefile_template = makefile_template.replace('vlasovsolver', 'vlasovsolver_amr')
    
with open('Makefile', 'w') as f:
    f.write(makefile_template)

# Finish with diagnostic output
display_settings()
