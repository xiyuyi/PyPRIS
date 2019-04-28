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
Upload: upload test_data, PyPRIS and .py files to hoffman2 under your home directory.

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

### Test the demo case on hoffman2:

Checkout a computation session with:

`qrsh`

Generate pypyrs configuration ticket with:

`python prep.py`

Run the PRIS recovery demo case with:

`python test.py`
