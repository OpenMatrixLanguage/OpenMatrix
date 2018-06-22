#!/bin/bash
    ## Script to launch the OpenMatrix 'omlconsole' on Linux platforms
	## Settings
	thedir=$(dirname $(readlink -m $BASH_SOURCE))
    export OML_ROOT=$thedir/..
    export PATH=$OML_ROOT/src/bin/linux64:$PATH
    export OML_THIRDPARTY=$OML_ROOT/../third_party
    export OML_HELP=$OML_ROOT/help/win/en/topics/reference/oml_language/
    # the user is advised to update the path as required on their system
    # add OpenMatrix binary directory 
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_ROOT/src/bin/linux64
    # Add fftw which needs to precede /lib64 directory from the system default
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/fftw/fftw-3.2.2/.libs:$LD_LIBRARY_PATH
    #add BLAS and LAPACK
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/lapack/lapack-3.7.1
    #add ANTLR
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/ANTLR/libantlr3c-3.4/.libs
    #add python3.5.2
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/python/python3.5
    #add HDF5 1.10.1
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/hdf/hdf5-1.8.16/shared/linux64/lib
    #add Matio 1.5.11
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/matio/matio-1.5.2_HDF/bin/linux64
    #add Sundials 3.1.0
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/sundials/sundials-3.1.0-install/lib
    #set pythonhome and pythonversion 
    export OML_PYTHONHOME=$OML_THIRDPARTY/python/python3.5
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OML_THIRDPARTY/python/python3.5/lib
	
	
echo Launch ATOM
# Update the line below (based on actual install location) to invoke ATOM 
atom


