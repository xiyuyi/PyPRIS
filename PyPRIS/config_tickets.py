import numpy as np
from PyPRIS import *

class BiplaneTicket:
    def __init__(self):
        self.plane1_path = None
        self.plane2_path = None
        self.psf_path = None
        "the name of this pris ticket"
        self.name = 'Demo'
        self.expansion = True
        "where to find the data files, for both blur and observation"
        self.datapath = './test_dataset_1'
        self.plane1_path = "{}/blur_plane1.tif".format(self.datapath)
        self.plane2_path = "{}/blur_plane7.tif".format(self.datapath)
        self.psf_path = "{}/psf.tif".format(self.datapath)
        self.psf_norm_factor = 80

        "specification of the psf stack"
        self.psfz0 = 106 # the count of the center plane, define it as z=0 in the psf coordinates system.
        self.plane1_dz = np.int8(-12)
        self.plane2_dz = np.int8(12)
        self.observer_edge_padding = True

        "configure the initial candidate pool of this pris ticket"
        self.init_candidates_intervals = list([1,5,5])
        self.init_ax0_range = list([-36, 36])
        self.init_ax1_range = list([5, 66])
        self.init_ax2_range = list([5, 66])

        "debug configurations"
        self.observer_debugger = False
        self.tobserver_edge_padding = True

        "output settings"
        self.output_file_header = 'Demo_biplane'
        self.output_path = './PyPRIS_Scratch'

        "linbreg configurations"
        self.linbreg_alpha = LinBreg("Demo")
        self.linbreg_alpha.debug = False
        self.linbreg_alpha.deep_debug = False
        self.linbreg_alpha.mu = 1000000000
        self.linbreg_alpha.alpha = 1e-8
        self.linbreg_alpha.maxit = 2000
        self.linbreg_alpha.it_check_rem = 1
        self.linbreg_alpha.debug_it_int = 100
        self.linbreg_alpha.kick.ints = 10
        self.linbreg_alpha.kick.flag = True
        self.linbreg_alpha.kick.thres = 1e-3
        self.linbreg_alpha.save_obj_int = 100
        self.linbreg_alpha.save = True
        self.linbreg_alpha.PyPRIS_iter = self.name
        self.linbreg_alpha.path_0 = self.output_path

        self.bg_scaling_coef = 1.5

        "others"
        self.PRIS_iter_end = 5

class SinglePlaneTicket:
    def __init__(self):
        self.plane1_path = None
        self.psf_path = None
        "the name of this pris ticket"
        self.name = 'Demo'
        self.expansion = True
        "where to find the data files, for both blur and observation"
        self.datapath = './test_dataset_2'
        self.plane1_path = "{}/im1.tif".format(self.datapath)
        self.psf_path = "{}/psf.tif".format(self.datapath)
        self.psf_norm_factor = 80

        "specification of the psf stack"
        self.psfz0 = 200 # the count of the center plane, define it as z=0 in the psf coordinates system.
        self.plane1_dz = np.int8(0)
        self.observer_edge_padding = True

        "configure the initial candidate pool of this pris ticket"
        self.init_candidates_intervals = list([1, 5, 5])
        self.init_ax0_range = list([-36, 36])
        self.init_ax1_range = list([5, 66])
        self.init_ax2_range = list([5, 66])

        "debug configurations"
        self.observer_debugger = False
        self.tobserver_edge_padding = True

        "output settings"
        self.output_file_header = 'Demo_single_plane'
        self.output_path = './PyPRIS_Scratch'

        "linbreg configurations"
        self.linbreg_alpha = LinBreg("Demo")
        self.linbreg_alpha.debug = False
        self.linbreg_alpha.deep_debug = False
        self.linbreg_alpha.mu = 1000000000
        self.linbreg_alpha.alpha = 1e-8
        self.linbreg_alpha.maxit = 2000
        self.linbreg_alpha.it_check_rem = 1
        self.linbreg_alpha.debug_it_int = 100
        self.linbreg_alpha.kick.ints = 10
        self.linbreg_alpha.kick.flag = True
        self.linbreg_alpha.kick.thres = 1e-3
        self.linbreg_alpha.save_obj_int = 100
        self.linbreg_alpha.save = True
        self.linbreg_alpha.PyPRIS_iter = self.name
        self.linbreg_alpha.path_0 = self.output_path
        self.bg_scaling_coef = 1.5

        "others"
        self.PRIS_iter_end = 5

class TwoChannelTicket:
    def __init__(self):
        "the name of this pris ticket"
        self.name = 'Demo'
        self.expansion = True

        "where to find the data files, for both blur and observation"
        # need to modify these default values with a prepared test dataset.
        self.datapath = './test_dataset_4'
        self.blur_path_channel_1 = "{}/BP+250_binAll.tif".format(self.datapath)
        self.blur_path_channel_2 = "{}/BP-250_binAll.tif".format(self.datapath)
        self.psf_path_channel_1 = "{}/psf_BP+250.tif".format(self.datapath)
        self.psf_path_channel_2 = "{}/psf_BP-250.tif".format(self.datapath)
        self.psf_norm_factor = 10000

        "specification of the psf stack"
        self.psfz0 = 76  # the count of the center plane, define it as z=0 in the psf coordinates system.
        self.observer_edge_padding = True

        "configure the initial candidate pool of this pris ticket"
        self.init_candidates_intervals = list([1, 5, 5])
        self.init_ax0_range = list([-60, 60])
        self.init_ax1_range = list([1, 64])
        self.init_ax2_range = list([1, 64])

        "debug configurations"
        self.observer_debugger = False
        self.tobserver_edge_padding = True

        "output settings"
        self.output_file_header = 'Demo_two_channel'
        self.output_path = './PyPRIS_Scratch'

        "linbreg configurations"
        self.linbreg_alpha = LinBreg("Demo")
        self.linbreg_alpha.debug = False
        self.linbreg_alpha.deep_debug = False
        self.linbreg_alpha.mu = 1000000000
        self.linbreg_alpha.alpha = 1e-8
        self.linbreg_alpha.maxit = 2000
        self.linbreg_alpha.it_check_rem = 1
        self.linbreg_alpha.debug_it_int = 100
        self.linbreg_alpha.kick.ints = 10
        self.linbreg_alpha.kick.flag = True
        self.linbreg_alpha.kick.thres = 1e-3
        self.linbreg_alpha.save_obj_int = 100
        self.linbreg_alpha.save = True
        self.linbreg_alpha.PyPRIS_iter = self.name
        self.linbreg_alpha.path_0 = self.output_path

        self.bg_scaling_coef = 1.5
        "others"
        self.PRIS_iter_end = 5
