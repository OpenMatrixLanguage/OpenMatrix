run_test () {
echo ==================== ====================  ==================== ====================
echo Run Core RegressionTests
cd $OML_ROOT/Tests/RegressionTests
echo
pwd
# update for travis
perl $OML_ROOT/Tests/regressOMLConsole.pl -travis
echo Completed core regression tests
echo ==================== ====================  ==================== ====================
display_files "*.log"
echo " "
echo Completed Core RegressionTests Log Display
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
    echo "Display RefLogs/$f file.."; 
    echo "--------------------"; 
    cat reflogs/$f; 
    echo "--------------------"; 
  done
}

run_test
