import sys
sys.path.append("/u/home/x/xiyuyi/bin")
sys.path.append("G:\\DH_localization\\PyPRIS")
from PyPRIS import *
import copy
import os

ticket = SinglePlaneTicket()

"the name of this pris ticket"
"ticket.name = 'Demo'"  # moved to loop

"where to find the data files, for both blur and observation"
paths = []
paths.append('G:\\DH_localization\\PyPRIS\\test_dataset_3')
ax0_ranges = []
ax0_ranges.append(list([-20, 20]))

ticket_folders = []
ticket_folders.append('PyPRIS_MT3D_Astig_fov5_bin_100')

for datapath, ax0_range, ticket_folder in zip(paths, ax0_ranges, ticket_folders):
    ticket.datapath = datapath
    ticket.observation_channel_n = 'SinglePlane'
    ticket.plane1_path = "{}/fov5.tif".format(ticket.datapath)
    ticket.psf_path = "{}/psf.tif".format(ticket.datapath)
    ticket.psf_norm_factor = 80

    "specification of the psf stack"
    ticket.psfz0 = 170  # the count of the center plane, define it as z=0 in the psf coordinates system.
    ticket.observer_edge_padding = True

    "configure the initial candidate pool of this pris ticket"
    ticket.init_candidates_intervals = list([1, 4, 4])
    ticket.init_ax0_range = ax0_range

    ticket.init_ax1_range = list([1, 61])
    ticket.init_ax2_range = list([1, 61])

    "debug configurations"
    ticket.observer_debugger = False
    ticket.tobserver_edge_padding = True

    "output settings"
    ticket.ticket_folder = ticket_folder

    "linbreg configurations"
    ticket.linbreg_alpha = LinBreg("X")

    # disable for hoffman2 tickets.

    # try:
    #     import matplotlib.pyplot as plt
    #     ticket.linbreg_alpha.debug = True
    #     ticket.linbreg_alpha.deep_debug = True
    # except RuntimeError:
    #     ticket.linbreg_alpha.debug = False
    #     ticket.linbreg_alpha.deep_debug = False

    ticket.linbreg_alpha.debug = True
    ticket.linbreg_alpha.deep_debug = False

    "ticket.linbreg_alpha.mu = 1000000000"  # move to loop
    "ticket.linbreg_alpha.alpha = 1e-11"  # move to loop
    ticket.linbreg_alpha.maxit = 40000
    ticket.linbreg_alpha.it_check_rem = 1
    ticket.linbreg_alpha.debug_it_int = 500
    ticket.linbreg_alpha.kick.ints = 1000
    ticket.linbreg_alpha.kick.eval_ints = 10
    ticket.linbreg_alpha.kick.flag = True
    ticket.linbreg_alpha.kick.thres = 0.01
    ticket.linbreg_alpha.save_obj_int = 1000
    ticket.linbreg_alpha.save = True
    ticket.expansion = True
    ticket.linbreg_alpha.PyPRIS_iter = 0
    "ticket.linbreg_alpha.PyPRIS_name = ticket.name"  # moved to loop
    ticket.linbreg_alpha.path_0 = '.'
    "ticket.bg_scaling_coef = 1.5 "  # moved to loop
    ticket.linbreg_alpha.stopping_loghistpercdelres_thres = -13

    "others"
    ticket.PRIS_iter_end = 6
    try:
        if not os.path.exists("../{}".format(ticket.ticket_folder)):
            os.mkdir("../{}".format(ticket.ticket_folder))
    except OSError:
        pass

    for bgSCF in list([1]):
        for mu in list([2e8]):
            for alpha in list([1e-10]):
                ticket_new = copy.deepcopy(ticket)
                ticket_new.name = "bgSCF" + str(bgSCF) + "_mu" + str("%1.1e" % mu) + "_alpha" + str("%1.1e" % alpha)
                ticket_new.bg_scaling_coef = copy.deepcopy(bgSCF)
                ticket_new.linbreg_alpha.PyPRIS_name = ticket_new.name
                ticket_new.linbreg_alpha.mu = mu
                ticket_new.linbreg_alpha.alpha = alpha
                try:
                    if not os.path.exists("../{}/{}".format(ticket_new.ticket_folder, ticket_new.name)):
                        os.mkdir("../{}/{}".format(ticket_new.ticket_folder, ticket_new.name))
                except OSError:
                    pass
                with open("../{}/{}/Go.pris_ticket".format(ticket_new.ticket_folder, ticket_new.name), "wb") as f:
                    pickle.dump(ticket_new, f, pickle.HIGHEST_PROTOCOL)
