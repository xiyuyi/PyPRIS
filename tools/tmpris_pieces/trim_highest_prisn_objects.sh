function tmpris.trim_highest_prisn_objects
{
  if [ $# -eq 0 ]
  then
    echo Syntax:
    echo tmpris.trim_highest_prisn_objects [option]
    echo This function should be called within a folder with a packet of tickets.
    echo each ticket folder should start with the str specified in the [option]
    echo [option] - str to store the pre-fix string of the folders that each contains an individual ticket.
    echo Example:
    echo tmpris.trim_highest_prisn_objects frame
    echo     This will detect all the folders in the present working direction with named frame*, and go into each folder,
    ecoh find the saved_objects folder in there, and find the objects stored in the saved_objects folder with the highest
    echo pris order, and delete them.
  fi

  for i in $(ls)
  do
  if [[ $i == frame*mu*bg* ]]
  then
      echo $i
      cd $i
          if [ -e ./saved_objects/done ]
          then
          :
          else
              if [ -e ./saved_objects ]
              then
                  cd saved_objects
                      # now find the highest prisn in the folder, and store it in the prisn var.
                      for prisn in 0 1 2 3 4 5 6
                      do
                          k=$(ls PyPRIS_pris"$prisn".file)
                          if [[ $k == PyPRIS_pris"$prisn".file ]]
                          then
                              currprisn=$prisn
                          fi
                      done
                      # now clean that highest prisn relevant objects.
                      rm *pris"$currprisn"*
                  cd ..
              fi
          fi
      cd ..
  fi
done
}