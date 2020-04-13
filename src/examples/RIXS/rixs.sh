#!/bin/bash
# Define feff command
pd=`pwd`
feffcmd="$pd/../../bin/feff"
feffdir="$pd/../../bin/Seq"
# 1. Locate the EDGE card and get the two edges. I am not being careful about comments after the card.
edges=( `grep -i '^[\ \t]*EDGE' feff.inp |sed -e's/EDGE//g' -e's/\r//g' -e's/\n//g'` )
for edge in ${edges[@]}
do
	# 3. Set up directories.
        echo "edge=$edge"
	mkdir -p ${edge}

	# 4. Make feff.inp for each edge
	if [ $edge == 'VAL' ] # Not being careful about case anywhere.
	then

		# Make VAL/feff.inp
		# Put COREHOLE NONE, RLPRINT, and EDGE into feff.inp
		echo 'COREHOLE NONE' > VAL/feff.inp
		echo 'RLPRINT' >> VAL/feff.inp
		echo "EDGE ${edges[0]}" >> VAL/feff.inp

		# Make sure that there is a XANES card. Grid doesn't matter, so just replace what is there.
		echo 'XANES 20' >> VAL/feff.inp

		# take out COREHOLE, RIXS, XES, EDGE, and HOLE cards. Use sed to change dos2unix just in case.
		#grep -v 'COREHOLE' feff.inp |grep -v 'RIXS' |grep -v 'XES' |grep -v 'EDGE' |grep -v 'HOLE'  >> VAL/feff.inp
		grep -v 'COREHOLE' feff.inp |grep -v 'RIXS' |grep -v 'XES' |grep -v 'EDGE' \
			|grep -v 'HOLE' |grep -v 'XANES' |sed -e's/\r//g' >> VAL/feff.inp
		cd VAL
                echo 'Calculating XANES for valence'
		$feffcmd >& VAL.log
		cd ..
		# Make XES/fef.inp
		mkdir -p XES		
		echo 'COREHOLE NONE' > XES/feff.inp
		echo "EDGE ${edges[0]}" >> XES/feff.inp
		
		# take out EGRID card and all of it's options, as well as XANES, EDGE, RIXS, and HOLE
		grep -v 'COREHOLE' feff.inp |grep -v 'RIXS' |grep -v 'XANES' |grep -v 'EDGE' \
			|grep -v 'HOLE' |grep -v 'EGRID' |grep -v '_grid' |sed -e's/\r//g' >> XES/feff.inp
		cd XES
		echo 'Calculating XES'
		$feffcmd >& XES.log
		cd ..
	else
		# put in EDGE, RLPRINT, and ICORE cards
		echo "EDGE $edge" > $edge/feff.inp
		echo "RLPRINT" >> $edge/feff.inp
		case ${edges[0]} in
		
		K) icore=1 
			;;
		L1) icore=2 
			;;
		L2) icore=3
			;;
		L3) icore=4
			;;
		M1) icore=5
			;;
		M2) icore=6
			;;
		M3) icore=7
			;;
		M4) icore=8
			;;
		M5) icore=10
			;;
		N1) icore=11
			;;
		N2) icore=12
			;;
		N3) icore=13
			;;
		N4) icore=14
			;;
		N5) icore=15
			;;
		*)  echo 'Unknown edge.'
		    exit
			;;
		esac
		echo "ICORE $icore" > $edge/feff.inp
		
		# Make $edge/feff.inp
		# Put COREHOLE NONE, RLPRINT, and EDGE into feff.inp
		echo 'COREHOLE RPA' >> $edge/feff.inp
		echo 'RLPRINT' >> $edge/feff.inp
		echo "EDGE ${edges[0]}" >> $edge/feff.inp

		# Make sure that there is a XANES card. Grid doesn't matter, so just replace what is there.
		echo 'XANES 20' >> $edge/feff.inp
		# take out XES, RIXS, HOLE, EDGE
		grep -v 'RIXS' feff.inp |grep -v 'XANES' |grep -v 'XES' |grep -v 'EDGE' \
			|grep -v '^\ *\t*HOLE' | grep -v 'EDGE' |sed -e's/[\n\r]//g' >> $edge/feff.inp
		cd $edge
		echo "Calculating $edge XANES"
		$feffcmd >& $edge.log
		cd ..
	fi
done

# 6. Run rixs module in main directory. 
echo 'Calculating RIXS'
$feffdir/rdinp >& rixs.log
$feffdir/atomic >> rixs.log
$feffdir/rixs >> rixs.log
