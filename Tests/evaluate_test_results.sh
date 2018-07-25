eval_regression_test() {

result=$(find $1 -maxdepth 1 -name '*.log' | wc -l)
result_string=" - OK"

if [ $result -gt 0 ]; then
  test_result=1;
  ls -alF $1/*.log
  result_string=" - THIS TEST HAS FAILED"
fi

echo Test result fails for $1 are $result $result_string

return $result
}
 
eval_toolbox_test() {

result=$(find toolboxes/$1/tests -maxdepth 1 -name '*.log' | wc -l)
result_string=" - OK"

if [ $result -gt 0 ]; then
  test_result=1;
  ls -alF toolboxes/$1/tests/*.log
  result_string=" - THIS TEST HAS FAILED"
fi

echo Test result fails for $1 are $result  $result_string

return $result
}

# setup default test_result.  If no tests have failed, this will be the return result
test_result=0

#Check core regression tests for the existance of any *.log files
eval_regression_test RegressionTests

#Check each toolbox test directory for the existance of any *.log files
eval_toolbox_test omlCAE
eval_toolbox_test omlCalculus
eval_toolbox_test omlDiffEq
eval_toolbox_test omlMathUtils
eval_toolbox_test omlMatio
eval_toolbox_test omlOptimization
eval_toolbox_test omlPolynom
#eval_toolbox_test omlpythonbridge
eval_toolbox_test omlSignals
eval_toolbox_test omlStatistics

echo test_result $test_result
if [ $test_result -ne 0 ]; then
  echo Note:  Failed tests fails the build
fi
exit $test_result
