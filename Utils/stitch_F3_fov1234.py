import sys
sys.path.append("G:/DH_localization/PyPRIS") # for local tests.
sys.path.append("~/Desktop/PyPRIS")
from PyPRIS import *
import os
print("")
print("")
print("")
# stitch 4 field of views together.
datapath_key_minimum_seed="/Users/xiyuyi/Desktop/PyPRIS/Documents/Data_key_minimum_seed"

for prisN in np.arange(0,8):
    exec("pris"+str(prisN)+"_names = list()")
    exec("pris"+str(prisN)+"_linbregs = list()")
    exec("pris"+str(prisN)+"_pyprises = list()")

# define the names of the linbreg objects for each pris iteration on fov1.
pris0_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris0_7001.file")
pris1_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris1_34520.file")
pris2_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris2_40000.file")
pris3_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris3_40000.file")
pris4_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris4_40000.file")
pris5_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris5_40000.file")
pris6_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris6_40000.file")
pris7_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris7_40000.file")

# define the names of the linbreg objects for each pris iteration on fov2.
pris0_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris0_4028.file")
pris1_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris1_40000.file")
pris2_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris2_40000.file")
pris3_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris3_40000.file")
pris4_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris4_40000.file")
pris5_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris5_40000.file")
pris6_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris6_40000.file")
pris7_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris7_8879.file")

# define the names of the linbreg objects for each pris iteration on fov3.
pris0_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris0_1696.file")
pris1_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris1_40000.file")
pris2_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris2_40000.file")
pris3_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris3_40000.file")
pris4_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris4_40000.file")
pris5_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris5_40000.file")
pris6_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris6_40000.file")
pris7_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris7_40000.file")

# define the names of the linbreg objects for each pris iteration on fov4.
pris0_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris0_1330.file")
pris1_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris1_25573.file")
pris2_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris2_40000.file")
pris3_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris3_40000.file")
pris4_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris4_40000.file")
pris5_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris5_40000.file")
pris6_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris6_40000.file")
pris7_names.append("PyPRIS_bgSCF2_mu1.0e+04_alpha1.0e-09_pris7_17227.file")


for Ind in np.arange(0, 4):
    pris0_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris0_names[Ind])
    pris1_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris1_names[Ind])
    pris2_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris2_names[Ind])
    pris3_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris3_names[Ind])
    pris4_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris4_names[Ind])
    pris5_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris5_names[Ind])
    pris6_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris6_names[Ind])
    pris7_linbregs.apend(datapath_key_minimum_seed + "/fov" + str(Ind + 1) + "/" + pris7_names[Ind])

