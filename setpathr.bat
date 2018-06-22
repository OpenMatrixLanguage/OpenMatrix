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
set OML_PYTHONHOME=c:/Python36
set OML_PYTHONVERSION=python36
set OML_PYTHON_NUMPYDIR=%OML_PYTHONHOME%/Lib/site-packages/numpy/core/include/numpy
@echo off

@echo on
set OML_THIRDPARTY=\oss\third_party

rem add lapack directory to path
set path=%OML_THIRDPARTY%\lapack\lapack-3.7.1-build\bin;%path%

rem add fftw directory to path
set path=%OML_THIRDPARTY%\fftw\fftw-3.2.2\fftw-3.2.2-libs\x64\Release;%path%

rem add matio directory to path
set path=%OML_THIRDPARTY%\matio\matio-1.5.11\visual_studio\x64\Release;%path%

rem add sundials directory to path
set path=%OML_THIRDPARTY%\sundials\sundials-3.1.0-install\lib;%path%

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
