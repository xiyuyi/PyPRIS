#!/u/home/x/xiyuyi/.conda/envs/PyPRIS_env/bin/python

print("")
print("")
print("--------------------------------------------------------")
print("")
print("")
print("")
print("PyPRIS: [add some description]")
print("Developers: Xiyu Yi, Xingjia Wang @ UCLA, 2019")
print("PI: Shimon Weiss")
print("[some more description]")

import sys
sys.path.append("/u/home/x/xiyuyi/bin")
from PyPRIS import *
print("")
print("")
print("")
print("Linbreg continue!")
path = "/u/scratch/x/xiyuyi/PyPRIS_tickets_set4/bgSCF1.5_mu1.0e+10_alpha1.0e-09/saved_objects"  # specifie datafile position
PyPRIS_name = "PyPRIS_bgSCF1.5_mu1.0e+10_alpha1.0e-09_pris0_160001"  # specify datafile name
PyPRIS_SensMx_name = "PyPRIS_bgSCF1.5_mu1.0e+10_alpha1.0e-09_pris0_SensingMx"  # specify datafile name
linbreg = loadCSSolver(path, PyPRIS_name, PyPRIS_SensMx_name)
linbreg.go()