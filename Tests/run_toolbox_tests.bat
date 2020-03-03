@echo off
echo Run toolbox tests

call :run_test omlCAE
call :run_test omlCalculus
call :run_test omlDiffEq
call :run_test omlGeometry
call :run_test omlMathUtils
call :run_test omlMatio
call :run_test omlOptimization
REM call :run_test omlpythonbridge
call :run_test omlPolynom
call :run_test omlSignals
call :run_test omlStatistics
goto :end

:run_test 
echo ==================== ====================  ==================== ====================
echo Run %1 toolbox tests
cd %OML_ROOT%/Tests/toolboxes/%1/tests
echo.
pwd
%OML_ROOT%/Tests/regressOMLConsole.pl
echo Completed %1 toolbox tests
echo ==================== ====================  ==================== ====================
if not exist *.log echo ALL %1 toolbox tests PASS! & echo. & goto :EOF
echo.
for %%a in ("*.log") do call :run_log_compare %%a
goto :EOF

:run_log_compare
echo Display %1 and reflogs\%1 files
echo -------------------- -------------------- -------------------- 
echo %1
echo -------------------- -------------------- -------------------- 
type %1
echo -------------------- -------------------- -------------------- 
echo.
echo reflogs\%1
echo -------------------- -------------------- -------------------- 
type reflogs\%1
echo -------------------- -------------------- -------------------- 
echo.
goto :EOF

:end
