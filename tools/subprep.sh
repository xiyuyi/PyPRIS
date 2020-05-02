#!/bin/bash
des="suball"
scp ./slurmhead ./suball
for i in $(ls)
do
  if [[ -d $i ]]
  then
  k="cd "$(pwd)"/"$i
  echo $k >> $des
  echo "chmod u+x pris" >> $des
  echo "srun -N1 -n1 pris & " >> $des
  echo "" >> $des
  fi

done
echo "wait" >> $des
echo "" >> $des
echo "echo '================ system output of the executable ends here ================'" >> $des
echo "echo 'finish time'" >> $des
echo "date" >> $des
