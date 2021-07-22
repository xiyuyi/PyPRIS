function tmpris.trim_highest_prisn_objects_here
{
  if [ $# -eq 0 ]
  then
    echo Syntax:
    echo tmpris.trim_highest_prisn_objects_here
    echo This function should be called within a ticket folder.
    echo Example:
    echo tmpris.trim_highest_prisn_objects_here go
    echo     This will detect all the folders in the present working direction with named frame*, and go into each folder,
    ecoh find the saved_objects folder in there, and find the objects stored in the saved_objects folder with the highest
    echo pris order, and delete them.
  fi

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

}