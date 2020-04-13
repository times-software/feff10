#!/bin/tcsh -f

set topdir = `pwd`
set inputfile = feff.inp
set report = $topdir/test-report

set dir1 = ~/science/feff-distro/JFEFF/feff90/examples/
set dir2 = ~/science/feff9/examples/

set list = `ls */*/$inputfile */*/*/$inputfile | sed -e "s/$inputfile//g" |grep -v REFERENCE`

rm $report
touch $report
echo "SUMMARY OF TEST RESULTS:" > $report
echo "========================" >> $report
echo "" >> $report
echo "Test                        Started            Finished  Crash   Avg.Deviation   Max.Deviation" >> $report
echo "----------------------------------------------------------------------------------------------" >> $report
foreach testcase ($list )

   echo "Starting test $testcase     `date`" 
   diff $dir1/$testcase/feff.inp $dir2/$testcase/feff.inp  

end

echo
echo
echo
cat $report
