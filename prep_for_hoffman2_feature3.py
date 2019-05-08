from PyPRIS import *
import copy
import os

ticket = BiplaneTicket()

"the name of this pris ticket"
"ticket.name = 'Demo'" # moved to loop

"where to find the data files, for both blur and observation"
ticket.datapath = '/u/scratch/x/xiyuyi/PyPRIS_data/S7_fov578_v4_set1_ox1165_oy1369'
ticket.plane1_path = "{}/DH_plane8.tif".format(ticket.datapath)
ticket.plane2_path = "{}/DH_plane14.tif".format(ticket.datapath)
ticket.psf_path = "{}/psf.tif".format(ticket.datapath)
ticket.psf_norm_factor = 80

"specification of the psf stack"
ticket.psfz0 = 106 # the count of the center plane, define it as z=0 in the psf coordinates system.
ticket.plane1_dz = np.int8(-12)
ticket.plane2_dz = np.int8(12)
ticket.observer_edge_padding = True

"configure the initial candidate pool of this pris ticket"
ticket.init_candidates_intervals = list([1,6,6])
ticket.init_ax0_range = list([-50, 80])
ticket.init_ax1_range = list([1, 81])
ticket.init_ax2_range = list([11, 81])

"debug configurations"
ticket.observer_debugger = False
ticket.tobserver_edge_padding = True

"output settings"
ticket.ticket_folder= 'PyPRIS_tickets_F3_set2_v4Set1'

"linbreg configurations"
ticket.linbreg_alpha = LinBreg("X")

# disabled for hoffman2 tickets.
#
# try:
#     import matplotlib.pyplot as plt
#     ticket.linbreg_alpha.debug = True
#     ticket.linbreg_alpha.deep_debug = False
# except RuntimeError:
#     ticket.linbreg_alpha.debug = False
#     ticket.linbreg_alpha.deep_debug = False

ticket.linbreg_alpha.debug = False
ticket.linbreg_alpha.deep_debug = False

"ticket.linbreg_alpha.mu = 1000000000" # move to loop
"ticket.linbreg_alpha.alpha = 1e-11" # move to loop
ticket.linbreg_alpha.maxit = 600000
ticket.linbreg_alpha.it_check_rem = 1
ticket.linbreg_alpha.debug_it_int = 100
ticket.linbreg_alpha.kick.ints = 10
ticket.linbreg_alpha.kick.flag = True
ticket.linbreg_alpha.kick.thres = 1e-3
ticket.linbreg_alpha.save_obj_int = 5000
ticket.linbreg_alpha.save = True
ticket.linbreg_alpha.PyPRIS_iter = 0
"ticket.linbreg_alpha.PyPRIS_name = ticket.name" # moved to loop
ticket.linbreg_alpha.path_0 = '.'
"ticket.bg_scaling_coef = 1.5 "  # moved to loop
ticket.linbreg_alpha.stopping_loghistpercdelres_thres = -15

"others"
ticket.PRIS_iter_end = 5
try:
    if not os.path.exists("../{}".format(ticket.ticket_folder)):
        os.mkdir("../{}".format(ticket.ticket_folder))
except OSError:
    pass

for bgSCF in list([1.5]):
    for mu in list([1e5, 1e6, 1e7]):
        for alpha in list([1e-10, 1e-11]):
            ticket_new = copy.deepcopy(ticket)
            ticket_new.name = "bgSCF"+str(bgSCF)+"_mu"+str("%1.1e"%mu)+"_alpha"+str("%1.1e"%alpha)
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
