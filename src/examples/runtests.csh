#!/bin/tcsh -f

set topdir = `pwd`
set inputfile = feff.inp
set feffprogram = $topdir/../bin/feff
set report = $topdir/test-report

unset analyzeonly
#set analyzeonly

set list = `ls */*/$inputfile */*/*/$inputfile | sed -e "s/$inputfile//g" |grep -v REFERENCE`
#set list = `ls EXAFS/*/$inputfile | sed -e "s/$inputfile//g" |grep -v REFERENCE`
#set list = "EXELFS/Cu/ "

if ($# ) then
   echo Taking examples list from input
   set list = $*
endif

rm $report
touch $report
echo "SUMMARY OF TEST RESULTS:" > $report
echo "========================" >> $report
echo "" >> $report
echo "Test                        Started            Finished  Crash   Avg.Deviation   Max.Deviation" >> $report
echo "----------------------------------------------------------------------------------------------" >> $report
foreach testcase ($list )

   # Run FEFF for the testcase
   echo "Starting test $testcase     `date`" 
   #echo -n "$testcase test started   `date '+%m/%d/%y %H:%M:%S'`" >> $report 
   printf "%-28s" $testcase >> $report
   printf "%-17s  " "`date '+%m/%d/%y %H:%M:%S'`"  >> $report 
   cd $topdir
   cd $testcase
   if (! $?analyzeonly ) $feffprogram
   echo Finished test $testcase     `date`
   #echo -n "   --  Finished  "  `date '+%H:%M:%S'` >> $report
   printf "%-8s"   `date '+%H:%M:%S'` >> $report

   # Did any modules crash?
   if ( -e ".feff.error" && ! -z ".feff.error" ) then
      echo "   THIS TEST CRASHED"
      echo "   Showing .feff.error :"
      cat ".feff.error"
      printf "%5s   "  'Yes' >> $report
   else
      echo No crash
      #echo "   --   No crash " >> $report
      printf "%5s   "  'No' >> $report
   endif

   # Try to quantify the difference in the spectrum compared to the standard.
   # Calculate   1/ N_energy  SUM (energy points)  [ |( mu(E) - mu_reference(e)) / mu_reference(E)|  ]
   # And x 100 to express as %
   if ( -e eels.dat ) then
      set newspectrum = `cat eels.dat | grep -v '#' | sed -e ' s/^ *// ; s/ *$// ; s/ \{1,\}/ /g' | cut -f2 -d' ' `
      set oldspectrum = `cat reference_eels.dat | grep -v '#' | sed -e ' s/^ *// ; s/ *$// ; s/ \{1,\}/ /g' | cut -f2 -d' ' `
   else if ( -e compton.dat ) then
      set newspectrum = `cat compton.dat | grep -v '#' | sed -e ' s/^ *// ; s/ *$// ; s/ \{1,\}/ /g' | cut -f2 -d' ' `
      set oldspectrum = `cat reference_compton.dat | grep -v '#' | sed -e ' s/^ *// ; s/ *$// ; s/ \{1,\}/ /g' | cut -f2 -d' ' `
   else
      set newspectrum = `cat xmu.dat | grep -v '#' | sed -e ' s/^ *// ; s/ *$// ; s/ \{1,\}/ /g' | cut -f4 -d' ' `
      set oldspectrum = `cat referencexmu.dat | grep -v '#' | sed -e ' s/^ *// ; s/ *$// ; s/ \{1,\}/ /g' | cut -f4 -d' ' `
   endif
   # Sometimes, defaults for energy mesh range etc. change, resulting in cut off energy mesh.  In that case, code below works.
   # However, if remaining energy points are not identical, rubbish is bound to happen - we don't attempt to fix that.
   set npointsnew = `echo $#newspectrum  | bc `
   set npointsold = `echo $#oldspectrum  | bc `
   if ( $npointsnew < $npointsold ) then
      set npoints = $npointsnew
   else
      set npoints = $npointsold
   endif
   #echo $npointsold $npointsnew $npoints
   set rdiff = 0.0
   set rdifftermmax = 0.0
   set i = 1 
   while (  $i <= $npoints )
      #Notice the extreme scale=30 below - it's pretty normal for regular eels to be 10^-15
      set munew=`echo $newspectrum[$i] | sed 's/E/\\*10\\^/' | sed 's/+//'` #because bc is fucking retarded
      set muold=`echo $oldspectrum[$i] | sed 's/E/\\*10\\^/' | sed 's/+//' `
      if (! ("$munew" == "$muold") )  then
         # ( ) around muold in denominator are necessary because bc doesn't understand scientific notation - it's seeing literal products
         set rdiffterm = `echo "scale=30;(( $munew - $muold) / ( $muold ) )" | bc`
         if ( `echo " $rdiffterm < 0 " | bc ` ) set rdiffterm = `echo "scale=10; - $rdiffterm "|bc` 
         #set rdiff = `echo "scale=30;$rdiff + (( $munew - $muold) / ( $muold ) )" | bc`
         set rdiff = `echo "scale=30;$rdiff + $rdiffterm " | bc`
         #echo "  " rdiffterm $rdiffterm rdiff $rdiff
         if ( `echo " $rdiffterm > $rdifftermmax " | bc ` ) set rdifftermmax = $rdiffterm
      endif
      #echo " $i $munew $muold " >> $topdir/checkmu
      @ i ++
   end
   set rdifftermmax = `echo "scale=2;$rdifftermmax * 100" | bc`
   set rdiff = `echo "scale=30;$rdiff / $npoints * 100.0" | bc`
   # Reduce the number of decimals to something legible:
   set rdiffprint = `echo  $rdiff | xargs printf "%8.3f"`
   set rdifftermmaxprint = `echo  $rdifftermmax | xargs printf "%8.2f"`
   echo "   Test is $rdiffprint % inaccurate in the spectrum.  (Max  $rdifftermmaxprint % off.)"
   #echo "   $rdiffprint % average deviation (max.  $rdifftermmaxprint % )" >> $report
   printf "%8.3f%%      %8.2f%%\n" $rdiff $rdifftermmax >> $report
end

echo
echo
echo
cat $report
