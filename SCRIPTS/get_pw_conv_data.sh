        for ((i=300; i<=900; i=i+50));do 
                #store final energies (only makes sense for runs
		#which follow sequentially:
                E=`tail -1 ${i}eV/OSZICAR` ; echo $i $E  >> PW_SUMMARY.dat
		F=`grep drift ${i}eV/OUTCAR` ; echo $i $F >> Forces_SUMMARY.dat
		#come back out of current dir
	
        done
