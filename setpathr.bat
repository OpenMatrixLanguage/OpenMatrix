@echo off
rem Set paths and enviroment variables to build OpenMatrix - release mode
rem run this script from the OML root directory
if not /%OML_ROOT%/==// call :error1 "OML_ROOT already defined!" & goto :end

@echo on
set OML_ROOT=%cd%
set path=%OML_ROOT%\src\bin\win64;%path%
set OML_ROOT=%OML_ROOT:\=/%
set OML_DEBUG=
rem set to local python installation directory 
rem set OML_PYTHONHOME=
rem set OML_PYTHONVERSION=
rem set OML_PYTHON_NUMPYDIR=
@echo off

@echo on
set OML_THIRDPARTY=\oss\third_party


rem add MKL directory for libiomp5md.dll
set path=%OML_THIRDPARTY%\intel\compilers_and_libraries_2019.5.281\windows\redist\intel64_win\compiler;%path%
rem add MKL directory 
set path=%OML_THIRDPARTY%\intel\compilers_and_libraries_2019.5.281\windows\redist\intel64_win\mkl;%path%

rem MKL related environment variables.
set KMP_AFFINITY=compact,1,0,granularity=fine

rem add fftw directory to path
set path=%OML_THIRDPARTY%\fftw\fftw-3.2.2\fftw-3.2.2-libs\x64\Release;%path%

rem add matio directory to path
rem DEL set path=%OML_THIRDPARTY%\matio\matio-1.5.12\visual_studio\x64\Release;%path%
set path=%OML_THIRDPARTY%\matio\matio-1.5.19\win64\bin;%path%

rem add hdf5 directory to path
set path=%OML_THIRDPARTY%\hdf\hdf5-1.12.0\win64\bin;%path%
rem add sundials directory to path

rem add qhull directory to path
set path=%OML_THIRDPARTY%\qhull\qhull-2015.2\bin;%path%

set OML_THIRDPARTY=%OML_THIRDPARTY:\=/%
set OML_HELP=%OML_ROOT%/help/win/en/topics/reference/oml_language/
@echo off
echo Display OML environment variables 
set OML
pause

rem set screen to black with white text
color 07
rem title %CD%
goto :end

:error1
echo %1 
goto :EOF

:end
