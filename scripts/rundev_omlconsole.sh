##!/bin/bash
if [[ ! -v OML_ROOT ]]; then
    export PATH=$PWD/src/bin/linux64:$PATH
    # the Linux demostration build documents the third party files in this directory.
    # Update this environmennt variable if third party files are built elsewhere.
    export OML_THIRDPARTY=~/oss/third_party
    export OML_ROOT=$PWD
    export OML_HELP=$PWD/help/win/en/topics/reference/oml_language/
    # the user is advised to update the path as required on their system
    # add OpenMatrix binary directory 
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src/bin/linux64
    # Add fftw directory which needed to precede an existing, older version on the development system
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/fftw/fftw-3.2.2/.libs:$LD_LIBRARY_PATH
    #add BLAS and LAPACK
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/lapack/lapack-3.7.1
    #add ANTLR
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/ANTLR/libantlr3c-3.4/.libs
	#add HDF5 1.10.1
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/hdf5/hdf5-1.10.1/hdf5/lib
    #add Matio 1.5.11
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/matio/matio-1.5.11/src/.libs
    #add Sundials 3.1.0
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/sundials/sundials-3.1.0-install/lib
    # the user is advised to update the path as required on their system
    # add Python shared library directory 
    #add python 3.5.2
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/python/python3.5/lib
	# update this path to the actual location of Python
	export OML_PYTHONHOME=/usr/python/python3.5
    export OML_PYTHONVERSION=python352

else
    echo "OML_ROOT is already defined: $OML_ROOT"
fi

echo Launch omlconsole
omlconsole -f scripts/loadtoolboxes.oml -continue $@



