run_test () {
echo ==================== ====================  ==================== ====================
echo Run $1 toolbox tests
cd $OML_ROOT/Tests/toolboxes/$1/tests
echo
pwd
# update for travis
perl $OML_ROOT/Tests/regressOMLConsole.pl -travis
echo Completed $1 toolbox tests
echo ==================== ====================  ==================== ====================
display_files "*.log"
echo " "
echo Completed $1 Log display
echo ==================== ====================  ==================== ====================
}

display_files() {
shopt -s nullglob
for f in $1; 
  do 
    echo " "; 
    echo "Display $f file.."; 
    echo "--------------------"; 
    cat $f; 
    echo "--------------------"; 
    echo " "; 
    echo "Display reflogs/$f file.."; 
    echo "--------------------"; 
    cat reflogs/$f; 
    echo "--------------------"; 
  done
}

run_test omlCAE
run_test omlCalculus
run_test omlDiffEq
run_test omlGeometry
run_test omlMathUtils
run_test omlMatio
run_test omlOptimization
run_test omlPolynom
run_test omlpythonbridge
run_test omlSignals
run_test omlStatistics

