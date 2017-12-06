!/bin/bash

DIR=/work/apm13/HydratePaper/New_Optimisations_28_4_2015/Optimisations_with_12d_PbElectrons_9_12_2015/PbI2/Convergence/PW
        #for ((i=50; i<=600; i=i+50));do
        for ((i=650; i<=750; i=i+50));do
       		echo $i
		mkdir $DIR/${i}eV
		cp $DIR/INCAR $DIR/KPOINTS $DIR/POSCAR $DIR/POTCAR $DIR/${i}eV/
		cp $DIR/VaspConv $DIR/${i}eV/VaspConv${i}
		cd $DIR/${i}eV
		echo $PWD
                #Read what ENCUT currently is:
                #j=$(grep ENCUT INCAR | cut -d "=" -f2)
                
		#change the ENCUT and save tmp copy:
                sed -i.tmp "s/.*ENCUT.*/ENCUT=${i}/g" INCAR
                
		#change job script for current dir
                #NB:Use '#' as seperator to overcome "/" escaping
                #   issues.
                sed -i.tmp "10s#.*DIR=.*#DIR=$PWD#g" VaspConv${i}
                
		#run VASP:
                qsub VaspConv${i}
                #store final energies (only makes sense for runs
		#which follow sequentially:
                #E=`tail -1 OSZICAR` ; echo $i $E  >> PW_SUMMARY.dat
        done
        echo "In $PWD"
done
