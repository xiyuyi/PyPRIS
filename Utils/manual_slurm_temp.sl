#!/bin/bash
##### These lines are for Slurm
#SBATCH -N 2
#SBATCH -p pbatch
#SBATCH -t 01:30:00
#SBATCH -A cancer
#SBATCH -J f13-1

##### These are shell commands
echo 'startig time'
date
id=$SLURM_JOBID
echo job id is $id
fpath=a
fname=a
date
echo '================ system output of the executable starts here ================'

for jobn in $(seq 10): # try to do ten jobs a node.  so request N nodes such taht 10*N is larger than the total job number.
do
  k="fpath=\$fpath$jobn"; eval $k
  k="fname=\$fname$jobn"; eval $k
  echo submitting $fpath/$fname
    cd $fpath
    chmod u+x $fname
    srun -N1 -n1 $fname &
  
done
  
 wait
  
echo Done
echo finish time 
date
