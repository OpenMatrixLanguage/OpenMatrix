@echo off
goto :argument_loop
:help_list
@echo off
rem help list and script documentation
@echo ---------- ---------- ---------- ---------- ---------- 
echo appveyor_test.bat  List of command line argument invoked options
echo mingw            - display contents of the mingw directory 
echo where            - run the windows where command for selected program files
echo end              - terminate this script without running the regression tests
echo core_test_perl   - run perl script regression tests
echo core_test_python - run python script regression tests
echo dir_x64_Release  - display exe/dll programs built in x64\Release directory 
echo check_x64_Release - check that all DLL/EXE are built in the x64\Release directory 
echo help             - display this help list
echo.
echo Note:  When OML_APPVEYOR is defined the core_test_perl regression tests are run by default.
echo.
rem
rem Summary of command line options:
rem appveyor_setup mingw where core_test_perl core_test_python end
rem 
goto :EOF

rem ---------- ---------- ---------- ---------- ---------- 
:argument_loop
if /%1/==// goto :break
if /%1/==/mingw/ call :appveyor_mingw & shift &goto :argument_loop
if /%1/==/where/ call :where_check & shift &goto :argument_loop
if /%1/==/core_test_perl/ call :core_test_perl & shift &goto :argument_loop
if /%1/==/core_test_python/ call :core_test_python & shift &goto :argument_loop
if /%1/==/appveyor_setup/ call :appveyor_setup & shift &goto :argument_loop
if /%1/==/dir_x64_Release/ call :dir_x64_Release & shift &goto :argument_loop
if /%1/==/check_x64_Release/ call :prog_x64_Release & shift &goto :argument_loop
if /%1/==/help/ call :help_list & shift &goto :argument_loop
if /%1/==/end/ goto :end
echo %1 - unrecognized command line argument
shift
goto :argument_loop

:break

@echo off
if defined OML_APPVEYOR call :appveyor_setup
if not defined OML_APPVEYOR echo There is no default option when OML_APPVEYOR is not defined. & echo Run script with 'help' argument to see command line options. & goto :end
call :core_test_perl
goto :end

:appveyor_setup
if exist appveyor_setup echo appveyor_setup already run & goto :EOF
@echo on 
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Run appveyor_setup
set appveyor_path=C:\mingw-w64\x86_64-7.2.0-posix-seh-rt_v5-rev1\mingw64\bin;
set appveyor_path=%OML_THIRDPARTY%\lapack\lapack-3.7.1-build\bin;%appveyor_path%
set appveyor_path=%OML_ROOT%\VS2015\OpenMatrix\x64\Release;%appveyor_path%
set path=%appveyor_path%;%path%
set appveyor_setup=TRUE
@echo off
goto :EOF 

:core_test_perl
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Run regression tests using perl script
cd %OML_ROOT%\Tests\RegressionTests
perl %OML_ROOT%\Tests\regressOMLConsole.pl
goto :EOF

:core_test_python
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Run regression tests using python script
cd %OML_ROOT%\Tests\RegressionTests
%OML_ROOT%\Tests\oml_run_tests.py -x omlconsole.exe
goto :EOF

:appveyor_mingw
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Display mingw directory on appveyor
dir C:\mingw-w64\x86_64-7.2.0-posix-seh-rt_v5-rev1\mingw64\bin
if errorlevel 1 echo errorlevel %errorlevel%
goto :EOF

:where_check
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Run the where command for key files
where omlconsole.exe
where oml.dll
where mathcore.dll
where libblas.dll
where liblapack.dll
goto :EOF

:dir_x64_Release
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Display programs built in x64\Release directory
dir c:\oss\OpenMatrix\VS2015\OpenMatrix\x64\Release\*.exe
dir c:\oss\OpenMatrix\VS2015\OpenMatrix\x64\Release\*.dll
goto :EOF

:check_x64_Release
@echo ---------- ---------- ---------- ---------- ---------- 
@echo Check that programs are built in the x64\Release directory 
set appveyor_found_count=0
set appveyor_notfound_count=0

for %%a in (omlconsole.exe oml.dll mathcore.dll hwcae.dll hwcalculus.dll hwdiffeq.dll hwmathutils.dll hwoptimization.dll hwpolynom.dll hwsignals.dll hwstatistics.dll   omlbatchplot.dll omlbatchplotcore.dll omlCAE.dll omlCalculus.dll omlDiffEq.dll omlMathUtils.dll omlMatio.dll omlOptimization.dll omlPolynom.dll omlpythonbridge.dll omlSignals.dll omlStatistics.dll) do call :prog_check %%a
echo OML FOUND   program count %appveyor_found_count%
echo OML MISSING program count %appveyor_notfound_count%
set appveyor_found_count
set appveyor_notfound_count

goto :EOF

:prog_check
if not exist %1 echo NOT FOUND %1 & set /a appveyor_notfound_count = appveyor_notfound_count + 1
if exist %1     echo FOUND %1  & set /a appveyor_found_count = appveyor_found_count + 1
rem c:\oss\OpenMatrix\VS2015\OpenMatrix\x64\Release\%1
goto :EOF


rem ---------- ---------- ---------- ---------- ---------- 
rem omlconsole.exe -e "disp('omlconsole.exe is running')"
rem omlconsole.exe -r abs1.oml
rem omlconsole.exe -r abs1.oml> abs1.log
rem fc abs1.log reflogs\abs1.log
rem \Oss\OpenMatrix\Tests\oml_run_tests.py -x omlconsole.exe
if exist abs1.log type abs1.log
if exist zeros17.log type zeros17.log
if exist abs1.log powershell format-hex abs1.log
powershell format-hex RefLogs\abs1.log
if exist zeros17.log powershell format-hex zeros17.log
powershell format-hex RefLogs\zeros17.log
git config --list
rem perl \Oss\OpenMatrix\Tests\regressOMLConsole.pl
rem type abs1.log
rem file compare text
rem fc abs1.log reflogs\abs1.log
rem file compare binary
rem fc /b abs1.log reflogs\abs1.log

:end