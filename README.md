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
load anaconda:

upload test_data, PyPRIS and .py files to hoffman2 under the home directory.

start an interactive session, load anaconda:

`qrsh`

`module load anaconda`

Install the enviroenemtn to your home directory with:
`conda env create -f PyPRIS_env.yml
