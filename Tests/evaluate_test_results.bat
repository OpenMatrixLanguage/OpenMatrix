@echo off
rem setup default test_result.  If no tests have failed, this will be the return result
set evaluate_result=0
set evaluate_verbose=true
if /%1/==/q/ set evaluate_verbose=
if /%1/==/quiet/ set evaluate_verbose=

rem Check core regression tests for the existance of any *.log files
call :eval_regression_test RegressionTests


rem Check each toolbox test directory for the existance of any *.log files
call :eval_toolbox_test omlCAE
call :eval_toolbox_test omlCalculus
call :eval_toolbox_test omlDiffEq
call :eval_toolbox_test omlGeometry
call :eval_toolbox_test omlMathUtils
call :eval_toolbox_test omlMatio
call :eval_toolbox_test omlOptimization
call :eval_toolbox_test omlPolynom
rem call :eval_toolbox_test omlpythonbridge
call :eval_toolbox_test omlSignals
call :eval_toolbox_test omlStatistics

echo evaluate_result %evaluate_result%
if /i %evaluate_result% NEQ 0 (
	echo.
    echo Note:  Failed tests fails the build
)

exit /b %evaluate_result%

goto :end

rem =======================================================
:filecount
rem @echo off
set filecountstr=%1
set filecount=0
for %%a in (%filecountstr%) do set /a filecount += 1
goto :EOF

rem =======================================================
:eval_regression_test
pushd %1
set filestr=%1
call :filecount *.log
if /i %filecount% EQU 0 (
    set result_str= - OK
)
if /i %filecount% GTR 0 (
	if defined evaluate_verbose (
		echo.
		dir *.log
		echo.
		for %%a in (*.log) do echo Failed %filestr%\%%a
		echo.
	)
    set result_str= - THIS TEST HAS FAILED
    set evaluate_result=1
)

echo Test result fails for %1 are %filecount% %result_str%
popd
goto :EOF

 
:eval_toolbox_test
pushd toolboxes\%1\tests
call :filecount *.log
if /i %filecount% EQU 0 (
    set result_str= - OK
)
if /i %filecount% GTR 0 (
	if defined evaluate_verbose (
		echo.
		dir *.log
		echo.		
	)
    set result_str= - THIS TEST HAS FAILED
    set evaluate_result=1
)

echo Test result fails for %1 are %filecount% %result_str%
popd
goto :EOF

:end