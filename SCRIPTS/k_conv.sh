#MAPI(Cubic) kpoints script.
#a=6.5, b=6.5, c=6.5
#so use i x i x i k points.
for (( i = 2; i <= 9; i=i+1 )) ; do  
echo "In ${PWD}"
echo "\n ${i}x${i}x${i}"
#make dir for run
mkdir k_${i}x${i}x${i}
#copy files to it
cp KPOINTS POSCAR POTCAR INCAR k_${i}x${i}x${i}
cp VaspConv k_${i}x${i}x${i}/MAPIK${i}x${i}x${i}
#change to created dir
cd k_${i}x${i}x${i}
# change KPOINTS (bihydrate is i x i x i):
sed -i.tmp "4s/.*[0-9].*/${i} ${i} ${i}/g" KPOINTS
#change job script for current dir
#NB:Use '#' as seperator to overcome "/" escaping
#   issues.
sed -i.tmp "10s#.*DIR=.*#DIR=$PWD#g" MAPIK${i}x${i}x${i}
# run VASP
qsub MAPIK${i}x${i}x${i}
#change back for next run 
cd ../
echo "In ${PWD}"
done
