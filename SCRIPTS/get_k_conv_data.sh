        for ((i=2; i<=8; i=i+2));do 
                a=3
b=2
h=$((a*i))
j=$((h / b))
		#store final energies (only makes sense for runs
		#which follow sequentially:
                E=`tail -1 k_${j}x${j}x${i}/OSZICAR` ; echo $j $j $i $E  >> K_SUMMARY.dat
		#come back out of current dir
	
        done
