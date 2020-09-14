@echo off
echo Current Directory is %cd%
pushd ..\..\..
set HW_ROOTDIR=%cd%
popd
pushd ..\RegressionTests
set HW_ROOTDIR
echo Test Directory is %cd%

perl regressVS.pl
popd

:end