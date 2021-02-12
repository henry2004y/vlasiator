#!/usr/bin/env python
# ------------------------------------------------------------------------------
# configure.py: Vlasiator configuration script in python.
#
# It uses the command line options and default settings to create customized
# versions of Makefile from the template file Makefile.in.
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

# -install
parser.add_argument('-install',
                    action='store_true',
                    default=False,
                    help='install all the dependencies into lib')

# --coord=[name]
parser.add_argument(
    '--coord',
    default='cartesian',
    choices=[
        'cartesian',
        'user'],
    help='select coordinate system')

# --nx=[value]
parser.add_argument('--nx',
                    type=int,
                    default=4,
                    help='set number of cells in the block x dimension')

# --ny=[value]
parser.add_argument('--ny',
                    type=int,
                    default=4,
                    help='set number of cells in the block y dimension')

# --nz=[value]
parser.add_argument('--nz',
                    type=int,
                    default=4,
                    help='set number of cells in the block z dimension')

# --velocityorder=[value]
parser.add_argument('--velocityorder',
                    type=int,
                    default=5,
                    help='set order of velocity space acceleration')

# --spatialorder=[value]
parser.add_argument('--spatialorder',
                    type=int,
                    default=3,
                    help='set order of spatial translation')

# --field=[name]
parser.add_argument('--field',
                    default='londrillo_delzanna',
                    help='select field solver')

# --fieldorder=[value]
parser.add_argument('--fieldorder',
                    type=int,
                    default=2,
                    help='select field solver order')

# -amr
parser.add_argument('-amr',
                    action='store_true',
                    default=False,
                    help='enable AMR in velocity space')

# -debug
parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable debug flags; override other compiler options')

# -debugsolver
parser.add_argument('-debugsolver',
                    action='store_true',
                    default=False,
                    help='enable debug flags for field solver')

# -debugionosphere
parser.add_argument('-debugionosphere',
                    action='store_true',
                    default=False,
                    help='enable debug flags for ionosphere module')

# -debugfloat
parser.add_argument('-debugfloat',
                    action='store_true',
                    default=False,
                    help='enable catching floating point exceptions')

# -float
parser.add_argument('-float',
                    action='store_true',
                    default=False,
                    help='enable single precision')

# -distfloat
parser.add_argument('-distfloat',
                    action='store_true',
                    default=True,
                    help='enable single precision for distribution function')

# -profile
parser.add_argument('-profile', 
                    dest='profile', 
                    action='store_true',
                    help='enable profiler')
parser.add_argument('-noprofile',
                    dest='profile',
                    action='store_false',
                    help='disable profiler')
parser.set_defaults(profile=True)

# -mpi
parser.add_argument('-mpi',
                    action='store_true',
                    default=True,
                    help='enable parallelization with MPI')

# -omp
parser.add_argument('-omp',
                    action='store_true',
                    default=False,
                    help='enable parallelization with OpenMP')

parser.add_argument('-nocheck',
                    action='store_true',
                    default=False,
                    help='disable linking library checking')


# --machine=[name]
parser.add_argument('--machine',
                    default='default',
                    help='choose machine name for specific setup')

# -papi
parser.add_argument('-papi',
                    action='store_true',
                    default=False,
                    help='enable Papi memory profiler')

# -jemalloc
parser.add_argument('-jemalloc',
                    action='store_true',
                    default=False,
                    help='enable JEMALLOC allocator')

# -silo
parser.add_argument('-silo',
                    action='store_true',
                    default=False,
                    help='enable silo format converter')

# The main choices for --cxx flag, using "ctype[-suffix]" formatting, where 
# "ctype" is the major family/suite/group of compilers and "suffix" may 
# represent variants of the compiler version and/or predefined sets of compiler
# options. The C++ compiler front ends are the main supported/documented options
# and are invoked on the command line, but the C front ends are also acceptable 
# selections and are mapped to the matching C++ front end:
# gcc -> g++, clang -> clang++, icc-> icpc
cxx_choices = [
    'g++',
    'g++-simd',
    'icpc',
    'icpc-debug',
    'icpc-phi',
    'cray',
    'clang++',
    'clang++-simd',
    'clang++-apple',
]


def c_to_cpp(arg):
    arg = arg.replace('gcc', 'g++', 1)
    arg = arg.replace('icc', 'icpc', 1)

    if arg == 'clang':
        arg = 'clang++'
    else:
        arg = arg.replace('clang-', 'clang++-', 1)
    return arg


# --cxx=[name]
parser.add_argument(
    '--cxx',
    default='g++',
    type=c_to_cpp,
    choices=cxx_choices,
    help='select C++ compiler and default set of flags')

# --cflag=[string]
parser.add_argument(
    '--cflag',
    default=None,
    help='additional string of flags to append to compiler/linker calls')

# Parse command-line inputs
args = vars(parser.parse_args())

# --- Step 2. Test for incompatible arguments ----------------------------

if args['install'] and args['machine'] is not 'default':
    raise SystemExit('### CONFIGURE ERROR: does not support fresh installation with a preset machine makefile')

# --- Step 3. Set Makefile options based on above argument

# Prepare dictionaries of substitutions to be made
makefile_options = {}
makefile_options['PREPROCESSOR_FLAGS'] = ''

if len(sys.argv) == 1:
    # Check existing Makefile if no argument is passed
    if not os.path.isfile("Makefile"):
        print("Vlasiator is not installed.")

# --cxx=[name] argument
if args['cxx'] == 'g++':
    makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -fopenmp -funroll-loops -ffast-math -std=c++17 -W -Wall -Wno-unused -mavx'
if args['cxx'] == 'g++-simd':
    # GCC version >= 4.9, for OpenMP 4.0; version >= 6.1 for OpenMP 4.5 support
    makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['COMPILER_FLAGS'] = (
        '-O3 -std=c++11 -fopenmp-simd -fwhole-program -flto -ffast-math '
        '-march=native -fprefetch-loop-arrays'
        # -march=skylake-avx512, skylake, core-avx2
        # -mprefer-vector-width=128  # available in gcc-8, but not gcc-7
        # -mtune=native, generic, broadwell
        # -mprefer-avx128
        # -m64 (default)
    )
if args['cxx'] == 'icpc':
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -inline-forceinline -qopenmp-simd -qopt-prefetch=4 '
      '-qoverride-limits'  # -qopt-report-phase=ipo (does nothing without -ipo)
    )
    # -qopt-zmm-usage=high'  # typically harms multi-core performance on Skylake Xeon
if args['cxx'] == 'icpc-debug':
    # Disable IPO, forced inlining, and fast math. Enable vectorization reporting.
    # Useful for testing symmetry, SIMD-enabled functions and loops with OpenMP 4.5
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -xhost -qopenmp-simd -fp-model precise -qopt-prefetch=4 '
      '-qopt-report=5 -qopt-report-phase=openmp,vec -g -qoverride-limits'
    )
if args['cxx'] == 'icpc-phi':
    # Cross-compile for Intel Xeon Phi x200 KNL series (unique AVX-512ER and AVX-512FP)
    # -xMIC-AVX512: generate AVX-512F, AVX-512CD, AVX-512ER and AVX-512FP
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xMIC-AVX512 -inline-forceinline -qopenmp-simd '
      '-qopt-prefetch=4 -qoverride-limits'
    )
if args['cxx'] == 'cray':
    makefile_options['COMPILER_COMMAND'] = 'CC'
    makefile_options['COMPILER_FLAGS'] = '-O3 -h std=c++14 -h aggress -h vector3 -hfp3'
    makefile_options['LINKER_FLAGS'] = '-hwp -hpl=obj/lib'
if args['cxx'] == 'clang++':
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++14'
if args['cxx'] == 'clang++-simd':
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++14 -fopenmp-simd'
if args['cxx'] == 'clang++-apple':
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++14'

# -float argument
if args['float']:
    precision = 'SP'
else:
    precision = 'DP'

makefile_options['PREPROCESSOR_FLAGS'] += ' -D'+precision

for key in ('vlsvdiff','vlsvextract','vlsv2silo'):
    makefile_options[key.upper()] = key+'_'+precision

# Distribution function precision
if args['distfloat']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DSPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4F_FALLBACK' # vector backend type
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4D_FALLBACK' # vector backend type

if args['profile']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DPROFILE'
 
if args['velocityorder'] == 2:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PLM'       
elif args['velocityorder'] == 3:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PPM'
elif args['velocityorder'] == 5:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PQM' 
else:
    raise SystemExit('### CONFIGURE ERROR: unknown semilag solver order for velocity space acceleration')

if args['spatialorder'] == 2:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PLM'       
elif args['spatialorder'] == 3:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PPM'
elif args['spatialorder'] == 5:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PQM' 
else:
    raise SystemExit('### CONFIGURE ERROR: unknown semilag solver order for spatial translation')

# Make the field solver first-order in space and time
if args['fieldorder'] == 1:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DFS_1ST_ORDER_SPACE'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DFS_1ST_ORDER_TIME'

# Turn on AMR
if args['amr']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DAMR'

# Add -DNDEBUG to turn debugging off.
if args['debug']:
    pass
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DNDEBUG'

# Debug for field solver
if args['debugsolver']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDEBUG_SOLVERS'

# Debug for ionosphere module
if args['debugionosphere']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDEBUG_IONOSPHERE'

# Catch floating point exceptions and stop execution
if args['debugfloat']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DCATCH_FPE'

if args['papi']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DPAPI_MEM'

if args['jemalloc']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DUSE_JEMALLOC -DJEMALLOC_NO_DEMANGLE'

# -debug argument
if args['debug']:
    # Completely replace the --cxx= sets of default compiler flags, disable optimization,
    # and emit debug symbols in the compiled binaries
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g'  # -Og
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] = '-O0 -h std=c++11'
    if args['cxx'] == 'bgxlc++':
        makefile_options['COMPILER_FLAGS'] = '-O0 -g -qlanglvl=extended0x'
    if args['cxx'] == 'icpc-phi':
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g -xMIC-AVX512'

# -mpi argument
if args['mpi']:
    if (args['cxx'] == 'g++' or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'icpc-phi' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_COMMAND'] = 'mpicxx'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -h mpi1'
else:
    raise SystemExit('### CONFIGURE ERROR: -mpi is required for compilation!')

# -omp argument
if args['omp']:
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd' or args['cxx'] == 'clang++'
            or args['cxx'] == 'clang++-simd'):
        makefile_options['COMPILER_FLAGS'] += ' -fopenmp'
    if (args['cxx'] == 'clang++-apple'):
        # Apple Clang disables the front end OpenMP driver interface; enable it via the
        # preprocessor. Must install LLVM's OpenMP runtime library libomp beforehand
        makefile_options['COMPILER_FLAGS'] += ' -Xpreprocessor -fopenmp'
    if args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi':
        makefile_options['COMPILER_FLAGS'] += ' -qopenmp'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -homp'
else:
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -hnoomp'
    if args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi':
        # suppressed messages:
        #   3180: pragma omp not recognized
        makefile_options['COMPILER_FLAGS'] += ' -diag-disable 3180'

# --cflag=[string]
if args['cflag'] is not None:
    makefile_options['COMPILER_FLAGS'] += ' '+args['cflag']

# Install dependencies
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
        errmsghead = """Boost not found: try to set the correct path, or """
        if 'Ubuntu' in os.uname()[3]:
            raise SystemExit(errmsghead+"install it manually as follows:"
            '\n sudo apt update\n sudo apt install libboost-all-dev')
        else:
            raise SystemExit(errmsghead+
            'search for how to install Boost!')

# --- Step 4. Create new files, finish up --------------------------------

# Read templates
with open('Makefile.in', 'r') as f:
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
print('Vlasiator has now been configured with the following options:')
print('  Machine:                    ' + (args['machine'] if args['machine'] else 'new'))
print('  Coordinate system:          ' + args['coord'])
print('  Floating-point precision:   ' + ('single' if args['float'] else 'double'))
print('  Distribution precision:     ' + ('single' if args['distfloat'] else 'double'))
print('  Block size:                 ' + str(args['nx']) + ' ' \
                                       + str(args['ny']) + ' ' \
                                       + str(args['nz']))
print('  MPI parallelism:            ' + ('ON' if args['mpi'] else 'OFF'))
print('  OpenMP parallelism:         ' + ('ON' if args['omp'] else 'OFF'))
print('  Order of field solver:      ' + str(args['fieldorder']))
print('  Order of semilag velocity:  ' + str(args['velocityorder']))
print('  Order of semilag spatial:   ' + str(args['spatialorder']))
print('  AMR:                        ' + ('ON' if args['amr'] else 'OFF'))
print('  Profiler:                   ' + ('ON' if args['profile'] else 'OFF'))
print('  Memory tracker:             ' + ('ON' if args['papi'] else 'OFF'))
print('  Debug flags:                ' + ('ON' if args['debug'] else 'OFF'))
print('  Compiler:                   ' + args['cxx'])
print('  Compilation command:        ' + makefile_options['COMPILER_COMMAND'] + ' '
      + makefile_options['COMPILER_FLAGS'])