import sys
sys.path.append("/u/home/x/xiyuyi/bin")
sys.path.append("../PyPRIS")
pypris_path = '../PyPRIS'
sys.path.append(pypris_path)
from PyPRIS import *
import copy
import os

ticket = TwoColor_CL_and_grating_Ticket()
# put field distortion parameters here. caused by CL:
# fitted parameters from pre-processing steps using pinhole array

"where to find the data files, for both blur and observation"
ticket.datapath = '/g/g92/yi10/PyPRIS/test_dataset_6'
ticket.blur_path_CL = "{}/view4_frame406_dif0CL.tif".format(ticket.datapath)
ticket.blur_path_diff1 = "{}/view4_frame406_dif1.tif".format(ticket.datapath)

ticket.psf_path_CL_Crimson = "{}/PSF_Crimson_Order0.tif".format(ticket.datapath)
ticket.distortion_Crimson_1_amp = 1
ticket.distortion_Crimson_1_shift = -0.421
ticket.distortion_Crimson_2_amp = 0.844
ticket.distortion_Crimson_2_shift = 17.484

ticket.psf_path_CL_Darkred = "{}/PSF_Darkred_Order0.tif".format(ticket.datapath)
ticket.distortion_Darkred_1_amp = 1
ticket.distortion_Darkred_1_shift = -0.421
ticket.distortion_Darkred_2_amp = 0.844
ticket.distortion_Darkred_2_shift = 17.484

ticket.psf_path_diff1_Crimson = "{}/PSF_Crimson_Order1.tif".format(ticket.datapath)
ticket.shift_Crimson_1 = 0.5
ticket.shift_Crimson_2 = 0

ticket.psf_path_diff1_Darkred = "{}/PSF_Darkred_Order1.tif".format(ticket.datapath)
ticket.shift_Darkred_1 = 2
ticket.shift_Darkred_2 = 0.2

ticket.psf_norm_factor = 80


"specification of the psf stack"
ticket.psfz0_Crimson = 101  # the count of the center plane, define it as z=0 in the psf coordinates system.
ticket.psfz0_Darkred = 101  # the count of the center plane, define it as z=0 in the psf coordinates system.
ticket.plane1_dz = np.int8(0)
ticket.observer_edge_padding = True

"configure the initial candidate pool of this pris ticket"
ticket.init_candidates_intervals = list([2, 4, 4])
ticket.init_ax0_range = list([-100, 100])

#ticket.init_ax1_range = list([1, 40])  # full range is 140  move to loop
#ticket.init_ax2_range = list([1, 40])  # full range is 212  move to loop


"debug configurations"
ticket.observer_debugger = False
ticket.observer_edge_padding = True

"output settings"
ticket.output_file_header = 'Demo_single_plane'
ticket.output_path = './PyPRIS_Scratch'

ticket.top_candidates = True
ticket.top_candidates_N = 100
"where to find the data files, for both blur and observation"

ticket.scratchpath = "/p/lscratchh/yi10/PRIS"
ticket_folders = []
ticket_folders.append('PyPRIS_V4F406_top100_tiles_50kIt')
ticket.ticket_folder = ticket_folders[0]


"linbreg configurations"
ticket.linbreg_alpha = LinBreg("X")
ticket.linbreg_alpha.debug = True
ticket.linbreg_alpha.deep_debug = False
"ticket.linbreg_alpha.mu = 1000000000"  # move to loop
"ticket.linbreg_alpha.alpha = 1e-11"  # move to loop


ticket.species_n = 2
ticket.linbreg_alpha.maxit = 1000000
ticket.linbreg_alpha.it_check_rem = 1
ticket.linbreg_alpha.debug_it_int = 5000
ticket.linbreg_alpha.kick.ints = 5000
ticket.linbreg_alpha.kick.eval_ints = 10
ticket.linbreg_alpha.kick.flag = True
ticket.linbreg_alpha.kick.thres = 0.01
ticket.linbreg_alpha.save_obj_int = 5000
ticket.linbreg_alpha.save = True
ticket.expansion = False
ticket.linbreg_alpha.PyPRIS_iter = 0
"ticket.linbreg_alpha.PyPRIS_name = ticket.name"  # moved to loop
ticket.linbreg_alpha.path_0 = '.'
"ticket.bg_scaling_coef = 1.5 "  # moved to loop
ticket.linbreg_alpha.stopping_loghistpercdelres_thres = -17

"others"
ticket.PRIS_iter_end = 7

try:
    if not os.path.exists("{}/{}".format(ticket.scratchpath, ticket.ticket_folder)):
        os.mkdir("{}/{}".format(ticket.scratchpath, ticket.ticket_folder))
except OSError:
    pass

ax1_ranges = [[1,45], [35,80], [70,115], [105,140]]
ax2_ranges = [[1,45], [35,80], [70,110], [105,140], [130, 175], [165,212]]

#ticket.init_ax1_range = list([1, 40])  # full range is 140  move to loop
#ticket.init_ax2_range = list([1, 40])  # full range is 212  move to loop
count = 0
for ticket.init_ax1_range in ax1_ranges:
    for ticket.init_ax2_range in ax2_ranges:
        for bgSCF in list([8]):
            for mu in list([1e2]):
                for alpha in list([1e-9]):
                    count += 1
                    ticket_new = copy.deepcopy(ticket)
                    ticket_new.name = "bgSCF" + str(bgSCF) + "_mu" + str("%1.1e" % mu) + "_alpha" + str("%1.1e" % alpha) + "_thres" + \
                                      str(ticket.linbreg_alpha.stopping_loghistpercdelres_thres) + "_c" + str(count).zfill(2) + \
                                      "_zrange" + str(ticket.init_ax0_range[0])+"to"+str(ticket.init_ax0_range[1]) + \
                                      "_yrange" + str(ticket.init_ax1_range[0]) + "to" + str(ticket.init_ax1_range[1]) + \
                                      "_xrange" + str(ticket.init_ax2_range[0]) + "to" + str(ticket.init_ax2_range[1])
                    ticket_new.bg_scaling_coef = copy.deepcopy(bgSCF)
                    ticket_new.linbreg_alpha.PyPRIS_name = ticket_new.name
                    ticket_new.linbreg_alpha.mu = mu
                    ticket_new.linbreg_alpha.alpha = alpha
                    try:
                        if not os.path.exists("{}/{}/{}".format(ticket_new.scratchpath, ticket_new.ticket_folder, ticket_new.name)):
                            os.mkdir("{}/{}/{}".format(ticket_new.scratchpath, ticket_new.ticket_folder, ticket_new.name))
                    except OSError:
                        pass
                    with open("{}/{}/{}/Go.pris_ticket".format(ticket_new.scratchpath, ticket_new.ticket_folder, ticket_new.name), "wb") as f:
                        pickle.dump(ticket_new, f, pickle.HIGHEST_PROTOCOL)
