#!/bin/bash
    ## Script to launch the OpenMatrix from within the ATOM Editor on Linux platforms
	## Settings
	thedir=$(dirname $(readlink -m $BASH_SOURCE))
    export OML_ROOT=$thedir/..
    export PATH=$OML_ROOT/src/bin/linux64:$PATH
    export OML_THIRDPARTY=$OML_ROOT/../third_party
    export OML_HELP=$OML_ROOT/help/win/en/topics/reference/oml_language/
    # the user is advised to update the path as required on their system
    # add OpenMatrix binary directory 
    export LD_LIBRARY_PATH=$OML_ROOT/src/bin/linux64:$LD_LIBRARY_PATH
    # Add fftw which needs to precede /lib64 directory from the system default
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/fftw/fftw-3.2.2/.libs:$LD_LIBRARY_PATH
    #add BLAS and LAPACK
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/lapack/lapack-3.7.1:$LD_LIBRARY_PATH
    #add ANTLR
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/ANTLR/libantlr3c-3.4/.libs:$LD_LIBRARY_PATH
    #add HDF5 1.10.1
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/hdf/hdf5-1.8.16/shared/linux64/lib:$LD_LIBRARY_PATH
    #add Matio 1.5.11
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/matio/matio-1.5.2_HDF/bin/linux64:$LD_LIBRARY_PATH
    #add Sundials 3.1.0
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/sundials/sundials-3.1.0-install/lib:$LD_LIBRARY_PATH
    #add zlib 1.2.8
    export LD_LIBRARY_PATH=$OML_THIRDPARTY/zlib-1.2.8/shared/lib/linux64:$LD_LIBRARY_PATH
    #set pythonhome and python library path 
    #export OML_PYTHONHOME=$OML_THIRDPARTY/python/python3.5
    #export LD_LIBRARY_PATH=$OML_THIRDPARTY/python/python3.5/lib:$LD_LIBRARY_PATH
    #export LD_LIBRARY_PATH=$OML_THIRDPARTY/python/python3.5:$LD_LIBRARY_PATH
	
echo Launch ATOM
# Update the line below (based on actual install location) to invoke ATOM 
atom


