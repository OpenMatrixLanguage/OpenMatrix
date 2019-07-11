@echo off
rem Set paths and environment variables to run OpenMatrix 'omlconsole' application

rem if environment variables are already defined skip to launch
if not /%OML_ROOT%/==// call :message1 "OML_ROOT already defined!" & goto :launch

rem @echo on
rem Run this script in the OML_ROOT directory to start omlconsole
set OML_ROOT=%cd%

set OML_DEBUG=
rem the Windows demonstration build documents the third party files in this directory.
rem Update this environment variable if third party files are built elsewhere.
set OML_THIRDPARTY=\oss\third_party
set OML_HELP=%OML_ROOT%/help/win/en/topics/reference/oml_language/

rem Set path to local installation of Python
set OML_PYTHONHOME=c:/Python36
set OML_PYTHONVERSION=python36

rem The user is advised to update the path as required on their system

rem add OpenMatrix executable directory to path
set path=%OML_ROOT%\src\bin\win64;%path%

rem Set paths to third parties used by OpenMatrix
rem add lapack directory to path
set path=%OML_THIRDPARTY%\lapack\lapack-3.7.1-build\bin;%path%
rem add fftw directory to path
set path=%OML_THIRDPARTY%\fftw\fftw-3.2.2\fftw-3.2.2-libs\x64\Release;%path%
rem add matio directory to path
set path=%OML_THIRDPARTY%\matio\matio-1.5.11\visual_studio\x64\Release;%path%
rem add sundials directory to path
set path=%OML_THIRDPARTY%\sundials\sundials-3.1.0-install\lib;%path%
rem add qhull directory to path
set path=%OML_THIRDPARTY%\qhull\qhull-2015.2\bin;%path%
rem add python directory to path
set path=%OML_PYTHONHOME%;%path%
rem add MINGW directory to path
set path=%OML_THIRDPARTY%\\mingw\x86_64-7.2.0-posix-seh-rt_v5-rev0\mingw\bin;%path%

@echo off

:launch
echo Launch omlconsole running script to load all toolboxes
omlconsole.exe -f "%OML_ROOT%"\scripts\loadtoolboxes.oml -continue
goto :end

:message1
echo %1 
goto :EOF

:end
