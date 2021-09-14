#!/bin/bash
#SBATCH --partition=batch  --ntasks=1 --cpus-per-task=48
#SBATCH --time=100:00:00
#SBATCH -x node7,node8
#SBATCH -e err
#SBATCH -o out
# record start date and host

 source /opt/intel/parallel_studio_xe_2020/bin/psxevars.sh intel64

export PATH=/opt/soft/vasp541:$PATH


if [ ! -f tag_0_finished ] ;then
  mpirun -np 48 /opt/soft/vasp/5.4.1/vasp_fix_z  1> runlog 2> errlog 
  if test $? -ne 0; then touch tag_failure_0; fi 
  touch tag_0_finished 
fi
