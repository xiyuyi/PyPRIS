function tmpris.countcalcpy
{

if [ $# -eq 9 ]
then
  slurmname=$2
  mem_per_task_G=$3
  mem_per_task=$((1024*$3))
  jobN=$4
  runtime=$5
  queue=$6
  bank=$7
  sequential_batch_N=$8
  folder_depth=$9
  if [[ $(hostname) == pascal* ]]
  then
    if [[ $queue == pdebug ]]
    then
    queue='pvis'
    fi
  fi
  unfinished=0
  total_tasks=0
  finished_tasks=0
  checkN=0
  for i in $(find . -maxdepth $folder_depth -type d)
  do
      echo $checkN
      checkN=`expr $checkN + 1`
          if [ -e $i/calcpy ]
          then
            total_tasks=`expr $total_tasks + 1`
            unfinished=`expr $unfinished + 1`
            k="task"$unfinished"_dir=\$(pwd)/\$i"
            eval $k
          fi
  done


  # find out total number of nodes
  if [[ $(hostname) == quartz* ]]
  then
    echo Youou are using quartz
    mem_per_node=128
  elif [[ $(hostname) == pascal* ]]
  then
    echo You are using pascal
    mem_per_node=256
  fi
  if [ $unfinished -ge 0 ]
  then
    # calculate tasks per job (round up formula), jobN is user input. one job is one slurm script.
    tasks_per_job=$((($unfinished+$jobN-1)/$jobN))

    # calculate tasks per job per batch (round up formula), sequential_batch_N is user input.
    tasks_per_job_per_batch=$((($tasks_per_job+$sequential_batch_N-1)/$sequential_batch_N))

    # calculate tasks per node based on memory requirement and memory availability (round-up formula)
    tasks_per_node=$(( ($mem_per_node) / $mem_per_task_G )) # this is fixed based on 2 user input value.
    # this is the parallel tasks per node.

    # calculate total number of nodes needed per job
    nodes_per_job=$((($tasks_per_job_per_batch+$tasks_per_node-1)/$tasks_per_node))

    # calculate the total number of tasks (from both sequential and parallel) needed to be calculated on one node.
    #    # use round up formula to get the maximum one.
    max_tasks_per_node=$((($tasks_per_job+$nodes_per_job-1)/$nodes_per_job))

    # calculate the actual batchN
    actual_batchN=$((($max_tasks_per_node+$tasks_per_node-1)/$tasks_per_node))

    # display status of task bundle
    echo You have a total of $total_tasks tasks under this directory and the subdirectories up to depth of $folder_depth
    echo $finished_tasks of them are finished
    echo $unfinished tasks to continue
    echo
    echo host: $(hostname)
    echo ---- memory per node: $mem_per_node
    echo ---- total jobs: $jobN  \(1 job corresponds to 1 slurm script\)
    echo ---- For each job: $nodes_per_job nodes will be requested, and will be used to run $tasks_per_job tasks in total.
    echo ---- For each node: up to $tasks_per_node tasks will be computed in the same instance.
    echo ---- $sequential_batch_N sequential batches are planned, actual batch number is up to $actual_batchN
    echo ---- each task reqiure $3 G memory.
    echo ---- Time to be requested for each job: $runtime
    echo ------- The total time is shared among all sequential tasks.
    echo ------- Please make sure the average time per batch is enough to finish one task.
    echo ----
    echo ---- tasks per job per sequential batch on each node is $tasks_per_job_per_batch
    echo ---- tasks per job on each node is $tasks_per_job
    echo ---- total number of batches is $sequential_batch_N
    echo
    echo
else
  echo all the tasks are finished.
  prepslurm=0
fi

  if [[ $1 == prepslurm ]]
  then
  tag=0
  batch_count=0
  # prepare for the slurm scripts
  echo "preparing slurm scripts... "
  for i in $(seq $jobN)
    do
      echo "prep job" $i ----- $slurmname"-"$i.sl
      sed -e "s/\[NODE NUMBER\]/$nodes_per_job/g" \
              -e "s/\[QUEUE TYPE\]/$queue/g" \
              -e "s/\[RUN TIME\]/$runtime/g" \
              -e "s/\[BANK\]/$bank/g"\
              -e "s/\[JOB NAME\]/$slurmname"-"$i/g" /g/g92/yi10/PyPRIS/tools/slurmhead > $slurmname"-"$i.sl

      for j in $(seq $tasks_per_job)
        do
          tag=`expr $tag + 1`
          batch_count=`expr $batch_count + 1`
          if [ $unfinished -ge $tag ]
            then
            k="task_dir=\$task"$tag"_dir"
            eval $k
            echo "echo submitting $task_dir/calcpy" >> $slurmname"-"$i.sl
            echo "  cd "$task_dir >> $slurmname"-"$i.sl
            echo "  chmod u+x calcpy" >> $slurmname"-"$i.sl
            echo "  srun -N1 -n1 --mem "$mem_per_task" calcpy &"   >> $slurmname"-"$i.sl
            echo "  " >> $slurmname"-"$i.sl
            if [ $batch_count -eq $tasks_per_job_per_batch ]
            then
              echo " wait" >> $slurmname"-"$i.sl
              echo "  " >> $slurmname"-"$i.sl
              batch_count=0
              fi
          fi
        done
        echo " wait " >>  $slurmname"-"$i.sl
        echo "echo '================ system output of the executable starts here ================'" >>  $slurmname"-"$i.sl
        echo "echo Done" >>  $slurmname"-"$i.sl
        echo "echo finish time " >>  $slurmname"-"$i.sl
        echo "date" >>  $slurmname"-"$i.sl
    done
  fi

else
  clear
  echo "tmpris.countcalcpy is very similar to tmpris.count, but it submit the calcpy script. It can do the following:"
  echo "    count how many calcpy scripts are there in the tree starting from the pwd"
  echo "    prepare slurmscripts"
  echo " you can choose the depth of folder-tree for the search."
  echo ""
  echo "options for input argument:"
  echo "    count: count the pack of tasks"
  echo "    prepslurm: to prepare a slurm script to submit all jobs"
  echo "usage:"
  echo " To count jobs:"
  echo "    tmpris.countcalcpy count <slurm-script-name> <mem-per-task(GB)> <jobN> <run time hh:mm:ss> <queue> <bank> <sequenctial_batch_N> <folder_depth>"
  echo " To prep slurm to submit jobs:"
  echo "    tmpris.countcalcpy prepslurm <slurm-script-name> <mem-per-task(GB)> <jobN> <run time hh:mm:ss> <queue> <bank> <sequenctial_batch_N> <folder_depth>"
  echo ""
  echo "Options"
  echo "<bank>:"
  echo "---- mmbp for LDRD-LW (project ended in 2020)"
  echo "---- cancer for spt-pris for ras"
  echo "---- dynimag for LDRD-ER"
  echo ""
  echo "<queue>:"
  echo "---- pbatch for calculation"
  echo "---- pdebug for debug"
  echo ""
  echo "<sequenctial_batch_N>:"
  echo "----  sequential batch number. It is sequential from batch to batch. tasks within the same batch are parallel"
fi


}