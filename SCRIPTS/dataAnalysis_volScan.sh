#!/bin/bash

home=$PWD
today=`date '+%d_%m_%Y'`;
for i in 0.90 0.95 1.00 1.05 1.10; do 
	# pad out the numbers 
	E=$(tail -1 Vol-${i}/OSZICAR)
	echo "${i} ${E}" > E_${i} 	
	P=$(grep 'pressure' Vol-${i}/OUTCAR)
	echo "${i} ${P}" > P_${i}
	cd $home
done
cat E_* > Energy.dat
rm E_*
cat P_* > Pressure.dat
rm P_*
#
