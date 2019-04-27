from PyPRIS import *
ticket = BiplaneConfig()

"the name of this pris ticket"
ticket.name = 'Demo'

"where to find the data files, for both blur and observation"
ticket.datapath = './test_data'
ticket.plane1_path = "{}/blur_plane1.tif".format(ticket.datapath)
ticket.plane2_path = "{}/blur_plane7.tif".format(ticket.datapath)
ticket.psf_path = "{}/psf.tif".format(ticket.datapath)
ticket.psf_norm_factor = 80

"specification of the psf stack"
ticket.psfz0 = 106 # the count of the center plane, define it as z=0 in the psf coordinates system.
ticket.plane1_dz = np.int8(-12)
ticket.plane2_dz = np.int8(12)
ticket.observer_edge_padding = True

"configure the initial candidate pool of this pris ticket"
ticket.init_candidates_intervals = list([1,5,5])
ticket.init_ax0_range = list([-36, 36])
ticket.init_ax1_range = list([5, 66])
ticket.init_ax2_range = list([5, 66])


"debug configurations"
ticket.observer_debugger = False
ticket.tobserver_edge_padding = True

"output settings"
ticket.output_file_header = 'Demo'
ticket.output_path = './DemoOutput'

"linbreg configurations"
ticket.linbreg_alpha = LinBreg("Demo")
ticket.linbreg_alpha.debug = True
ticket.linbreg_alpha.deep_debug = False
ticket.linbreg_alpha.mu = 1000000000
ticket.linbreg_alpha.alpha = 1e-8
ticket.linbreg_alpha.maxit = 200
ticket.linbreg_alpha.it_check_rem = 1
ticket.linbreg_alpha.debug_it_int = 100
ticket.linbreg_alpha.kick.ints = 10
ticket.linbreg_alpha.kick.flag = True
ticket.linbreg_alpha.kick.thres = 1e-3
ticket.linbreg_alpha.save_obj_int = 100
ticket.linbreg_alpha.save = True
ticket.bg_scaling_coef = 1.5

with open("{}.pris_ticket".format(ticket.name), "wb") as f:
    pickle.dump(ticket, f, pickle.HIGHEST_PROTOCOL)
