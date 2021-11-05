#!/bin/bash
#SBATCH --partition=batch  --ntasks=1 --cpus-per-task=48
#SBATCH --time=1000:00:00
#SBATCH -e err
#SBATCH -o out
# record start date and host
echo "date      --> `date`"
echo "hostname  --> `hostname`"
echo "total cpu --> $SLURM_NTASKS"
echo "name      -->  $SLURM_JOB_NAME"
echo "dir       -->  $SLURM_SUBMIT_DIR"
echo "nodelist  -->  $SLURM_JOB_NODELIST"

#  load the module, you may have to change this part
ulimit -s unlimited
module load vasp/5.4.1

# assume elc_energy or elc_stress folder in current directory
path=`pwd`
approach=elc_energy
cd $path/$approach

list1=`ls -d Def_*`
cd $path/$approach/Def_1
list2=`ls -d Def_1_*`
cd $path

for i in $list1
do
    echo --------- Def_$i -----------
    for j in $list2
    do
    # use the follow code to avoid duplicate submissions
       if [ ! -f tag_finished ] ;then
         # you may have to change the path of vasp command 
         mpirun -np 48 /opt/soft/vasp/5.4.1/vasp_fix_z  1> runlog 2> errlog 
         if test $? -ne 0; then touch tag_failure; echo Def_${i}_$j fail; fi 
         touch tag_finished 
         echo Def_${i}_$j ok
       fi
    done
done
  
