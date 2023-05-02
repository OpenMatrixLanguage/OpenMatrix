@echo off
rem The user is advised to update paths as required on their system.
rem Here we are updating paths to following the example in the BUILD_WINDOWS document.
rem Run this script from the OpenMatrix root directory.

rem if environment variables are already defined skip to launch
if not /%OML_ROOT%/==// call :message1 "OML_ROOT already defined!" & goto :launch

rem add OpenMatrix executable directory to path
set OML_ROOT=%cd%
set path=%OML_ROOT%\src\bin\win64;%path%
set OML_THIRDPARTY=%OML_ROOT%/../Continuous_Integration_Windows/third_party
set OML_HELP=%OML_ROOT%/help/win/en/topics/reference/oml_language/
set OML_DEBUG=

rem Set path to local installation of Python
rem set OML_PYTHONHOME=
rem set OML_PYTHONVERSION=
rem set OML_PYTHON_NUMPYDIR=

rem MKL related environment variables.
set KMP_AFFINITY=compact,1,0,granularity=fine

rem add MKL directory for libiomp5md.dll
set path=%OML_THIRDPARTY%\intel\compilers_and_libraries_2019.5.281\windows\redist\intel64_win\compiler;%path%
rem add MKL directory 
set path=%OML_THIRDPARTY%\intel\compilers_and_libraries_2019.5.281\windows\redist\intel64_win\mkl;%path%

rem add fftw directory to path
set path=%OML_THIRDPARTY%\fftw\fftw-3.2.2\fftw-3.2.2-libs\x64\Release;%path%

rem add matio directory to path
set path=%OML_THIRDPARTY%\matio\matio-1.5.19\win64\bin;%path%

rem add hdf5 directory to path
set path=%OML_THIRDPARTY%\hdf\hdf5-1.12.0\win64\bin;%path%

rem add qhull directory to path
set path=%OML_THIRDPARTY%\qhull\qhull-2015.2\bin;%path%


rem add MINGW directory to path
rem set path=%OML_THIRDPARTY%\\mingw\x86_64-7.2.0-posix-seh-rt_v5-rev0\mingw\bin;%path%

@echo off

:launch
echo Launch omlconsole running script to load all toolboxes
omlconsole.exe -f "%OML_ROOT%"\scripts\loadtoolboxes.oml -continue
goto :end

:message1
echo %1 
goto :EOF

:end
