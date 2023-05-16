#!/bin/bash

# The user is advised to update paths as required on their system.
# Here we are updating paths to following the example in BUILD_LINUX document.
# Run this script from the OpenMatrix root directory.
export PATH=$PWD/src/bin/linux64:$PATH
export OML_THIRDPARTY=~/ombuild/Continuous_Integration_Linux_OS66/third_party
export OML_ROOT=$PWD
export OML_HELP=$PWD/help/win/en/topics/reference/oml_language/

#Set python environment variables
#export OML_PYTHONHOME=
#export OML_PYTHONVERSION=
#export OML_PYTHON_NUMPYDIR=

# Save the system default
export LD_LIBRARY_PATH_SAVE=$LD_LIBRARY_PATH

export KMP_AFFINITY=compact,1,0,granularity=fine

# add directory for libiomp5md.dll
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin
export PATH=$PATH:$OML_THIRDPARTY/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin
# add MKL directory 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin
export PATH=$PATH:$OML_THIRDPARTY/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin

# Add fftw which needs to precede /lib64 directory from the system default
export LD_LIBRARY_PATH=$OML_THIRDPARTY/fftw/fftw-3.2.2/.libs:$LD_LIBRARY_PATH

#add OpenMatrix binary directory 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src/bin/linux64

#add ANTLR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/ANTLR/libantlr3c-3.4/.libs

#add Matio
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/matio/matio-1.5.19/linux64/lib64

#add HDF5
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/hdf/hdf5-1.12.0/linux64/lib

#add qhull 2015.2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/qhull/qhull-2015.2/lib

echo Display OML environment variables:
set | grep OML

echo -e "\033]0;<Configured for $OML_ROOT environment>\007"

