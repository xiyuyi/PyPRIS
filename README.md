# PyPRIS
This is a python package for PRIS -- a progressive refinement method for sparse recovery.

Authors: Xiyu Yi, Xingjia Wang @ UCLA

PI: Shimon Weiss

## Environment Setup

PyPRIS requires installation of Anaconda.

Before running PyPRIS, create the environment by running the following code in Anaconda:

`conda env create -f PyPRIS_env.yml`

Then activate the environment with:

`conda activate PyPRIS_env`

Before running notebook files, set up the environment of Ipython kernel as well with the following code:

`ipython kernel install --user --name=PyPRIS_env`

Finally, after starting jupyter notebook, switch kernel to current environment by clicking "Kernel -> Change kernel -> PyPRIS_env".

## On Hoffman2:
Upload: upload the entire test_data and PyPRIS folders, *PyPRIS_env.yml* and *\*.py* files to hoffman2 under your home directory.

#### configure environment for PyPRIS (PyPRIS_env) on Hoffman2 from terminal (recommended terminals include MobaXTerm and Putty) 
After log-in from the terminal, start an interactive session with: `qrsh` (this step may take a few minutes)

Load anaconda with:

`module load anaconda`

Install the PyPRIS_env environment to your home directory with:

`conda env create -f PyPRIS_env.yml` (This step may take about 10 minutes)

Activate the environment with:

`source activate PyPRIS_env`

To log-out from the PyPRIS environment, type: 

`conda deactivate`

### Computation on Hoffman2
* Upload the observation files (e.g. *blur_plane1.tif*, *blur_plane7.tif and *psf.tif* in test_data) to Hoffman2. Keep both files under the same root folder. 
* Prepare tickets by updating and executing *prep_for_hoffman2.py* locally (recommended.)
Examples include *prep_for_hoffman2_Feature2.py* and *prep_for_hoffman2_feature3.py*
  * Specify file names and file paths for each plane.
  * Specify ticket.ticket_folder
  * Specify index range for x,y,z dimensions (the range need to match your fov size)
  * Set a wide range for bgCSF, mu, and alpha values, in the form of absolute values.
  * Upload the entire ticket folders to Hoffman2.
  * Copy pris and pris.Ini files into each ticket folder.
  * Turn on or off **memory profiling options** in *pris* file 
      * Line 133: if a previously saved PyPRIS object already exists 
      * Line 178: for new calculations
      * Top 10 memory usages by line for each PyPRIS iteration will be saved in the saved_objects folder under the name *PyPRIS_*_mem.file*
  * To ensure the *pris* file, the ticket files and the associated file folders are readable, writable and executable, use `chmod -R u+rwx * `
  * **Perform single job test!** This is extremley important. Test with a single job first, as if you fail a big batch of jobs, your priority will be reduced in the queue.
      * If any error arises, open an interactive session to debug on the server by activating an interactive Python session. To initiate such session from the terminal on Hoffman2 (recommended to use PuTTY), type `Python`
      * Fix bugs until a linbreg recovery can be successfully executed in the interactive session.
      * **Perform single job test before moving on!**
      * Test with `python pris` from the terminal.
  * Submit batch jobs.
  * Remember that there are complications associated with files from Github that was pushed from a Windows system!
      * Temporary solution: creat the same file on the linux system by copy-pasting the content in to a new file with the same file name (instead of copy-pasting the file).
