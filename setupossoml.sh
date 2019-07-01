#!/bin/bash

export PATH=$PWD/src/bin/linux64:$PATH
export OML_THIRDPARTY=~/oss/third_party_os66
export OML_ROOT=$PWD
export OML_HELP=$PWD/help/win/en/topics/reference/oml_language/
export OML_PYTHONHOME=/usr
export OML_PYTHONVERSION=python3.6
export OML_PYTHON_NUMPYDIR=$OML_PYTHONHOME/lib64/$OML_PYTHONVERSION/site-packages/numpy/core/include/numpy/

# Save the system default
export LD_LIBRARY_PATH_SAVE=$LD_LIBRARY_PATH

# Add fftw which needs to precede /lib64 directory from the system default
export LD_LIBRARY_PATH=$OML_THIRDPARTY/fftw/fftw-3.2.2/.libs:$LD_LIBRARY_PATH

#add OpenMatrix binary directory 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src/bin/linux64

#add BLAS and LAPACK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/lapack/lapack-3.7.1

#add ANTLR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/ANTLR/libantlr3c-3.4/.libs

#add Matio 1.5.11
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/matio/matio-1.5.12/src/.libs

#add Sundials 3.1.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/sundials/sundials-3.1.0-install/lib

#add qhull 2015.2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/qhull/qhull-2015.2/lib

echo Display OML environment variables:
set | grep OML

echo -e "\033]0;<Configured for $OML_ROOT environment>\007"

