@echo off
echo Run tests

call :run_core_test
goto :end

:run_core_test 
echo ==================== ====================  ==================== ====================
echo Run core regression tests
cd RegressionTests
echo.
pwd
@echo on
perl %OML_ROOT%\Tests\regressOMLConsole.pl
@echo off
echo Completed core tests
echo ==================== ====================  ==================== ====================
if not exist *.log echo ALL core tests PASS! & echo. & goto :EOF
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
