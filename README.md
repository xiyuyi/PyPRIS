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
Upload: upload test_data, PyPRIS, PyPRIS_env.yml and \*.py files to hoffman2 under your home directory.

#### configure environment for PyPRIS (PyPRIS_env) on Hoffman2 from terminal (recommended terminals include MobaXTerm and Putty) 
After log-in from the terminal, start an interactive session with: `qrsh` (this step may take a few minutes)

Load anaconda with:

`module load anaconda`

Install the PyPRIS_env enviroenemt to your home directory with:

`conda env create -f PyPRIS_env.yml` (This step may take about 10 minutes)

Activate the environment with:

`source activate PyPRIS_env`

To log-out from the PyPRIS environment, type: 

`conda deactivate`

### Computation on Hoffman2
* upload the observation files (for example: blur_plane1.tif, blur_plane7.tif and psf.tif in test_data) to hoffman2. keep both files under the same file folder 
* prepare tickets by updating and executing prep_for_hoffman2.py locally (recommended.)
Examples include prep_for_hoffman2_Feature2.py and prep_for_hoffman2_feature3.py
  * specify file names and file paths for each plane.
  * specify ticket.ticket_folder
  * specify index range for x,y,z dimensions (the range need to match your fov size)
  * set a wide range for bgCSF, mu, and alpha values, in the form of absolute values.
  * upload to hoffman2.
  * put in the pris and pris.Ini files into each ticket folder.
  * To ensure the pris and the ticket files are readable and executable, use `chmod -R u+rwx * `
  * **Perform Single Job test!** extremley important. test with a single job first. Because if you fail a big batch of jobs, your priority will be reduced in the queue.
      * if any error arises, open an interactive session to debug on the server by activating an interactive python session.
          * to initiate an interactive python session from the terminal on hoffman2 (recommended to use Putty), type `Python`
      * fix bugs until a linbreg recovery can be successfully executed in the interactive session.
      * **perform single JOb test before moving on!**
      * test with `python pris` from the terminal.
  * submit batch jobs
