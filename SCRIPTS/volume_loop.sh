#!/bin/bash

home=$PWD
for i in 0.925 0.95 0.975 1.00 1.025 1.05 1.075; do 
	# pad out the numbers 
	mkdir Vol-$i
	cp vasp_single_point.qsub POSCAR POTCAR KPOINTS INCAR Vol-$i
	cd Vol-$i
	sed -i "2s/.*[0-9].*/${i}/g" POSCAR
	qsub vasp_single_point.qsub
	echo $i
	cd $home
done

#
