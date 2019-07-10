@echo off
rem Set paths and environment variables to run OpenMatrix
rem from within the ATOM Editor

rem if environment variables are already defined skip to launch
if not /%OML_ROOT%/==// call :message1 "OML_ROOT already defined!" & goto :launch

rem set this to start omlconsole from OML_ROOT
rem set OML_ROOT=%cd%

rem set this to start omlconsole from any folder
set OML_ROOT=%~dp0\..

set OML_DEBUG=
set OML_THIRDPARTY=%OML_ROOT%\..\third_party
set OML_HELP=%OML_ROOT%/help/win/en/topics/reference/oml_language/

rem Set path to local installation of Python
rem set OML_PYTHONHOME=%OML_THIRDPARTY%\python\python3.5

rem The user is advised to update the path as required on their system

rem add OpenMatrix executable directory to path
set path=%OML_ROOT%\src\bin\win64;%path%

rem Set paths to third parties used by OpenMatrix
rem add lapack directory to path
set path=%OML_THIRDPARTY%\lapack\lapack-3.7.1-build\bin;%path%
rem add fftw directory to path
set path=%OML_THIRDPARTY%\fftw\fftw-3.2.2\fftw-3.2.2-libs\x64\Release;%path%
rem add matio directory to path
set path=%OML_THIRDPARTY%\matio\matio-1.5.2_HDF\bin\win64;%path%
set path=%OML_THIRDPARTY%\hdf\hdf5-1.8.16\shared\win64\bin;%path%
rem add sundials directory to path
set path=%OML_THIRDPARTY%\sundials\sundials-3.1.0-install\lib;%path%
rem add qhull directory to path
set path=%OML_THIRDPARTY%\qhull\qhull-2015.2\bin;%path%
rem add python directory to path
rem set path=%OML_THIRDPARTY%\python\python3.5;%path%
rem add MINGW directory to path
set path=%OML_THIRDPARTY%\mingw\x86_64-7.2.0-posix-seh-rt_v5-rev0\mingw\bin;%path%
rem zlib 1.2.8
set path=%OML_THIRDPARTY%\zlib-1.2.8\shared\lib\win64;%path%
rem system_dlls
set path=%OML_THIRDPARTY%\system_dlls\vc14\win64;%path%


@echo off

:launch
echo Launch ATOM
rem Update the line below (based on actual install location) to invoke ATOM 
C:\Users\UserName\AppData\Local\atom\atom.exe

goto :end

:message1
echo %1 
goto :EOF

:end
