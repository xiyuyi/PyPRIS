import sys
sys.path.append("G:/DH_localization/PyPRIS") # for local tests.
from PyPRIS import *
import os
print("")
print("")
print("")
import warnings
warnings.filterwarnings("ignore")

# load the original ticket
ticket = get_ticket('./Go.pris_ticket')
ticket.datapath = 'G:/DH_localization/Experiments_Phase2_afterSPIE/April2019_BoulderTrip_analysis_choose_first_FOV_initialTest/S7_fov578_v4_set2_ox1155_oy1392'
ticket.plane1_path = "{}/DH_plane8.tif".format(ticket.datapath)
ticket.plane2_path = "{}/DH_plane14.tif".format(ticket.datapath)
ticket.psf_path = "{}/psf.tif".format(ticket.datapath)
ticket.linbreg_alpha.maxit = 40000
ticket.linbreg_alpha.save_obj_int = 200
ticket.linbreg_alpha.debug = True
ticket.linbreg_alpha.stopping_loghistpercdelres_thres = -11
ticket.expansion = True
with open("./Go_local.pris_ticket", "wb") as f:
    pickle.dump(ticket, f, pickle.HIGHEST_PROTOCOL)