#!/bin/bash

home=$PWD
today=`date '+%d_%m_%Y'`;
for i in 0.925 0.95 0.975 1.00 1.025 1.05 1.075; do 
	# pad out the numbers 
	E=$(tail -1 Vol-${i}/OSZICAR)
	V=$(grep -m 1 'volume of cell' Vol-${i}/OUTCAR | grep -o '[0-9]*\.[0-9]*')
	echo "${V} ${E} ${i}" > E_${i} 	
	P=$(grep 'pressure' Vol-${i}/OUTCAR)
	echo "${V} ${P}" > P_${i}
	cd $home
done
cat E_* > Energy.dat
rm E_*
cat P_* > Pressure.dat
rm P_*
#
