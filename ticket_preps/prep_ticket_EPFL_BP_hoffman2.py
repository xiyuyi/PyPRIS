import sys
sys.path.append("/u/home/x/xiyuyi/bin")
sys.path.append("G:\\DH_localization\\PyPRIS")

from PyPRIS import *
import copy
import os

ticket = TwoChannelTicket()

"where to find the data files, for both blur and observation"
paths = []
paths.append('/u/scratch/x/xiyuyi/EPFL_BP/test_dataset_4')
ax0_ranges = []
ax0_ranges.append(list([-60, 60]))
ticket_folders = []
ticket_folders.append('PyPRIS_EPFL_BP_binAll_hoffman2')

for datapath, ax0_range, ticket_folder in zip(paths, ax0_ranges, ticket_folders):
    "the name of this pris ticket"
    ticket.name = 'Demo'
    ticket.expansion = True

    "where to find the data files, for both blur and observation"
    # need to modify these default values with a prepared test dataset.
    ticket.datapath = datapath
    ticket.blur_path_channel_1 = "{}/BP+250_binAll.tif".format(ticket.datapath)
    ticket.blur_path_channel_2 = "{}/BP-250_binAll.tif".format(ticket.datapath)
    ticket.psf_path_channel_1 = "{}/psf_BP+250.tif".format(ticket.datapath)
    ticket.psf_path_channel_2 = "{}/psf_BP-250.tif".format(ticket.datapath)
    ticket.psf_norm_factor = 10000

    "specification of the psf stack"
    ticket.psfz0 = 76  # the count of the center plane, define it as z=0 in the psf coordinates system.
    ticket.observer_edge_padding = True

    "configure the initial candidate pool of this pris ticket"
    ticket.init_candidates_intervals = list([1, 3, 3])
    ticket.init_ax0_range = list([-75, 75])
    ticket.init_ax1_range = list([1, 64])
    ticket.init_ax2_range = list([1, 64])

    "debug configurations"
    ticket.observer_debugger = False
    ticket.tobserver_edge_padding = True

    "output settings"
    ticket.ticket_folder = ticket_folder

    "linbreg configurations"
    ticket.linbreg_alpha = LinBreg("X")
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
    ticket.expansion = False
    ticket.linbreg_alpha.PyPRIS_iter = 0
    "ticket.linbreg_alpha.PyPRIS_name = ticket.name"  # moved to loop
    ticket.linbreg_alpha.path_0 = '.'
    "ticket.bg_scaling_coef = 1.5 "  # moved to loop
    ticket.linbreg_alpha.stopping_loghistpercdelres_thres = -13

    "others"
    ticket.PRIS_iter_end = 8

    try:
        if not os.path.exists("../{}".format(ticket.ticket_folder)):
            os.mkdir("../{}".format(ticket.ticket_folder))
    except OSError:
        pass

    for bgSCF in list([2, 5]):
        for mu in list([5e6, 1e7, 5e7, 1e8]):
            for alpha in list([1e-4, 1e-5, 1e-6]):
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