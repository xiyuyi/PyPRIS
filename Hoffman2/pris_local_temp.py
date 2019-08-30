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

sys.path.append("/u/home/a/alejandr/bin")
sys.path.append("/u/home/x/xiyuyi/bin")
sys.path.append("G:\\DH_localization\\PyPRIS")
from PyPRIS import *
import os

print("")
print("")
print("")
import warnings

warnings.filterwarnings("ignore")

# load the original ticket
ticket = get_ticket('./Go.pris_ticket')

# prep the order according to the ticket [PSF, and overall blur]
if isinstance(ticket, SinglePlaneTicket):
    "load psf and blur for Single Plane Ticket"
    psf = io.imread(ticket.psf_path)
    blur = io.imread(ticket.plane1_path)

if isinstance(ticket, BiplaneTicket):
    "load psf and blur for sequential single channel biplane Ticket"
    psf = io.imread(ticket.psf_path)
    blur1 = io.imread(ticket.plane1_path)
    blur2 = io.imread(ticket.plane2_path)
    blur = np.concatenate([blur1, blur2], axis=0)

if isinstance(ticket, TwoChannelTicket):
    "load psf and blur for simultaneous two channel observation."
    psf_channel_1 = io.imread(ticket.psf_path_channel_1)
    psf_channel_2 = io.imread(ticket.psf_path_channel_2)
    blur_channel_1 = io.imread(ticket.blur_path_channel_1)
    blur_channel_2 = io.imread(ticket.blur_path_channel_2)
    blur = np.concatenate([blur_channel_1, blur_channel_2], axis=0)

if isinstance(ticket, SingleCLTicket):
    "load psf and blur for Single Plane Ticket"
    psf = io.imread(ticket.psf_path)
    blur = io.imread(ticket.blur_path)

if isinstance(ticket, CLnShiftTicket):
    "load psf and blur for ticket of CL channel and shift channel"
    psf_CL_in = io.imread(ticket.psf_path_CL)
    blur_CL_in = io.imread(ticket.blur_path_CL)
    psf_shift_in = io.imread(ticket.psf_path_shift)
    blur_shift_in = io.imread(ticket.blur_path_shift)
    blur = np.concatenate([blur_CL_in, blur_shift_in], axis=0)

# define linbreg_ori
linbreg_ori = ticket.linbreg_alpha
linbreg_ori.obs_dim0 = blur.shape[0]
linbreg_ori.obs_dim1 = blur.shape[1]
linbreg_ori.PyPRIS_name = ticket.name
linbreg_ori.stopping_loghistpercdelres_thres = -28
linbreg_ori.maxit = 140000
linbreg_ori.debug_it_int = 1000
# prepare a observer for pypris. [from one candidate, to one observation]
# summon an observer from the ObserveStation (ObserveStation class) (Cheers~)
observer = ObserveStation()
if isinstance(ticket, SinglePlaneTicket):
    "configure this observer with single plane (monoplane) observation specifics"
    observer.observe_monoplane_prep(psf, single_image_size=blur.shape,
                                    psfz0=ticket.psfz0,
                                    observer_debugger=ticket.observer_debugger,
                                    observer_edge_padding=ticket.observer_edge_padding)

if isinstance(ticket, BiplaneTicket):
    "configure this observer with biplane observation specifics"
    observer.observe_biplane_prep(psf, single_image_size=blur1.shape,
                                  deltaz_plane1=ticket.plane1_dz,
                                  deltaz_plane2=ticket.plane2_dz,
                                  psfz0=ticket.psfz0,
                                  observer_debugger=ticket.observer_debugger,
                                  observer_edge_padding=ticket.observer_edge_padding)

if isinstance(ticket, TwoChannelTicket):
    "configure this observer with two channel observation specifics"
    "simultaneous observation with two channels"
    observer.observe_two_channel_prep(psf1=psf_channel_1,
                                      size1=blur_channel_1.shape,
                                      psf2=psf_channel_2,
                                      size2=blur_channel_2.shape,
                                      psfz0=ticket.psfz0,
                                      observer_debugger=ticket.observer_debugger,
                                      observer_edge_padding=ticket.observer_edge_padding)

if isinstance(ticket, SingleCLTicket):
    "configure this observer with single observation that has a cylindrical lens"
    observer.observe_with_CL_prep(psf,
                                  single_iamge_size=blur.shape,
                                  psfz0=ticket.psfz0,
                                  dist_1_amplitd=ticket.distortion_1_amp,
                                  dist_1_shift=ticket.distortion_1_shift,
                                  dist_2_amplitd=ticket.distortion_2_amp,
                                  dist_2_shift=ticket.distortion_2_shift,
                                  observer_debugger=ticket.observer_debugger,
                                  observer_edge_padding=ticket.observer_edge_padding)

if isinstance(ticket, CLnShiftTicket):
    "configure this observer with both CL and shift channels"
    observer.observe_with_CL_and_grating_prep(psf_CL=psf_CL_in,
                                              imsize_CL=blur_CL_in.shape,
                                              psfz0_CL=ticket.psfz0_CL,
                                              dist_1_amplitd=ticket.distortion_1_amp,
                                              dist_1_shift=ticket.distortion_1_shift,
                                              dist_2_amplitd=ticket.distortion_2_amp,
                                              dist_2_shift=ticket.distortion_2_shift,
                                              psf_shift=psf_shift_in,
                                              imsize_shift=blur_shift_in.shape,
                                              psfz0_shift=ticket.psfz0_shift,
                                              shift_1=ticket.shift_channel_shift1,
                                              shift_2=ticket.shift_channel_shift2,
                                              observer_debugger=ticket.observer_debugger,
                                              observer_edge_padding=ticket.observer_edge_padding)

# Get the latest saved linbreg filename
# put in the file folder where you saved all the linbreg objects.
# path = 'C:\\Users\\wxjpp\\Desktop\\PyPRIS-master\\saved_objects'
path = './saved_objects'
mxnames = []
objIter = []

# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for file in f:
        if file.endswith('SensingMx.file'):
            mxnames.append(file[0:-5])

if len(mxnames) is not 0:
    print("PyPRIS Continue!")
    print("")
    mxName = mxnames[-1]
    for r, d, f in os.walk(path):
        for file in f:
            if file.startswith(mxName[0:-9]):
                if not file.endswith('SensingMx.file'):
                    objIter.append(int(file[len(mxName[0:-9]):-5]))
    #
    #
    objName = mxName[0:-9] + str(max(objIter))
    #
    print("Loading LinBreg object from: ")
    print(path)
    print(mxName)
    print(objName)
    #
    # Load the previous saved LinBreg object
    with open('{}/{}.file'.format(path, objName), "rb") as f:
        linbreg = pickle.load(f)  # the loaded object is a LinBreg object
    #
    with open('{}/{}.file'.format(path, mxName), "rb") as s:
        linbreg.A = joblib.load(s)
    #
    pypris = loadPyPRIS(path, 'PyPRIS_pris{}'.format(linbreg.PyPRIS_iter[4:]))
    # set a check mark for pypris:
    pypris.set_check_mark()
    # put the configured observer to the pypris to perform the tasks of 'observe_biplane'
    # when a pypris object observes, it is through the observer and it is
    # done with the observer's observe_biplane skill (skill = method).
    if isinstance(ticket, SinglePlaneTicket):
        "use the monoplane observation skill of the observer for pypris object"
        pypris.observe = observer.observe_monoplane

    if isinstance(ticket, BiplaneTicket):
        "use the biplane observation skill of the observer for pypris object"
        pypris.observe = observer.observe_biplane

    if isinstance(ticket, TwoChannelTicket):
        "use the two-channel observation skill of the observer for pypris object"
        pypris.observe = observer.observe_two_channel

    if isinstance(ticket, SingleCLTicket):
        "use the observe with single cylindrical lens observation skill of the observer for pypris object"
        pypris.observe = observer.observe_with_CL

    if isinstance(ticket, CLnShiftTicket):
        "use the observe with the combination of cylindrical lens channel observation skill and the shift channel for pypris object"
        pypris.observe = observer.observe_with_CL_and_grating

    Iter = int(linbreg.PyPRIS_iter[4:])
    print("---------------- PRIS refinement #" + str(Iter) + " ------------------")
    print("")
    linbreg.stopping_loghistpercdelres_thres = -28
    linbreg.maxit = 140000
    linbreg.debug_it_int = 1000
    linbreg.go()
else:
    print("PyPRIS Begin!")
    print("")
    # condition a main linbreg object that carries the basic configurations.
    # for each PRIS iteration, a new linbreg object will be cloned from this one for its own task.
    # the purpose here is just to simplify the parameter configuration process.
    linbreg_ori = ticket.linbreg_alpha
    linbreg_ori.obs_dim0 = blur.shape[0]
    linbreg_ori.obs_dim1 = blur.shape[1]
    linbreg_ori.PyPRIS_name = ticket.name
    #
    # get a PyPRIS object:
    pypris = PyPRIS()
    pypris.observation = blur.ravel()
    pypris.current_candidates_intervals = ticket.init_candidates_intervals  # initialize the first intervals of neighboring candidate voxels.
    #
    # decide candidates for the first PRIS iteration. PRIS #0:
    range_ind0 = np.arange(ticket.init_ax0_range[0], ticket.init_ax0_range[1], pypris.current_candidates_intervals[
        0])  # seems like for DH psf, use finnest intervals for the thickness dimension is better.
    range_ind1 = np.arange(ticket.init_ax1_range[0], ticket.init_ax1_range[1], pypris.current_candidates_intervals[1])
    range_ind2 = np.arange(ticket.init_ax2_range[0], ticket.init_ax2_range[1], pypris.current_candidates_intervals[2])
    #
    # initialize the first pool of candidates with the current candidate intervals.
    pypris.current_candidates = list()  # current pool of candidates.
    for i0 in range_ind0:
        for i1 in range_ind1:
            for i2 in range_ind2:
                pypris.current_candidates.append([i0, i1, i2])

    # set a check mark for pypris:
    #
    pypris.set_check_mark()
    pypris.save()
    # put the configured observer to the pypris to perform the tasks of 'observe_biplane'
    # when a pypris object observes, it is through the observer and it is
    # done with the observer's observe_biplane skill (skill = method).
    if isinstance(ticket, SinglePlaneTicket):
        "use the monoplane observation skill of the observer for pypris object"
        pypris.observe = observer.observe_monoplane

    if isinstance(ticket, BiplaneTicket):
        "use the biplane observation skill of the observer for pypris object"
        pypris.observe = observer.observe_biplane

    if isinstance(ticket, TwoChannelTicket):
        "use the two-channel observation skill of the observer for pypris object"
        pypris.observe = observer.observe_two_channel

    if isinstance(ticket, SingleCLTicket):
        "use the observe with single cylindrical lens observation skill of the observer for pypris object"
        pypris.observe = observer.observe_with_CL

    if isinstance(ticket, CLnShiftTicket):
        "use the observe with the combination of cylindrical lens channel observation skill and the shift channel for pypris object"
        pypris.observe = observer.observe_with_CL_and_grating

    # construct sensing matrix,
    pypris.generate_sensing_mx()
    # prepare the inner sparse recovery
    pypris.current_A = pypris.current_A / ticket.psf_norm_factor
    c = np.dot(pypris.observation.ravel(), pypris.current_A)
    c1 = np.min(c[0:-2])
    c2 = np.max(c[0:-2])
    bgv = c1 + (c2 - c1) * ticket.bg_scaling_coef
    pypris.current_A[:, -1] = bgv / np.sum(
        pypris.observation.ravel())  # the last column in A corresponds to the background component.

    linbreg = copy.deepcopy(linbreg_ori)
    linbreg.candidate_coords = pypris.current_candidates
    linbreg.candidate_intervals = pypris.current_candidates_intervals
    linbreg.A = pypris.current_A
    linbreg.b = pypris.observation.ravel()
    linbreg.PyPRIS_iter = "pris0"
    print("")
    print("---------------- PRIS refinement #0 ------------------")
    print("")
    # recover
    linbreg.get_ready()
    linbreg.go()
    Iter = 0

# now continue the pris iterations
# construct sensing matrix
for PRIS_iter in np.arange(Iter + 1, 9):
    pypris.prep_for_new_refinement()
    pypris.expansion = ticket.expansion
    pypris.refine_candidates(linbreg)
    pypris.current_PRIS_ItN += 1
    pypris.save()
    # put the configured observer to the pypris to perform the tasks of 'observe_biplane'
    # when a pypris object observes, it is through the observer and it is
    # done with the observer's observe_biplane skill (skill = method). or obsreve monoplane skill
    if isinstance(ticket, SinglePlaneTicket):
        "use the monoplane observation skill of the observer for pypris object"
        pypris.observe = observer.observe_monoplane

    if isinstance(ticket, BiplaneTicket):
        "use the biplane observation skill of the observer for pypris object"
        pypris.observe = observer.observe_biplane

    if isinstance(ticket, TwoChannelTicket):
        "use the two-channel observation skill of the observer for pypris object"
        pypris.observe = observer.observe_two_channel

    if isinstance(ticket, SingleCLTicket):
        "use the observe with single cylindrical lens observation skill of the observer for pypris object"
        pypris.observe = observer.observe_with_CL

    if isinstance(ticket, CLnShiftTicket):
        "use the observe with the combination of cylindrical lens channel observation skill and the shift channel for pypris object"
        pypris.observe = observer.observe_with_CL_and_grating

    pypris.generate_sensing_mx()
    # prepare the inner sparse recovery
    pypris.current_A = pypris.current_A / ticket.psf_norm_factor
    c = np.dot(pypris.observation.ravel(), pypris.current_A)
    c1 = np.min(c[0:-2])
    c2 = np.max(c[0:-2])
    bgv = c1 + (c2 - c1) * ticket.bg_scaling_coef
    pypris.current_A[:, -1] = bgv / np.sum(
        pypris.observation.ravel())  # the last column in A corresponds to the background component.
    linbreg = copy.deepcopy(linbreg_ori)
    linbreg.candidate_coords = pypris.current_candidates
    linbreg.candidate_intervals = pypris.current_candidates_intervals
    linbreg.A = pypris.current_A
    linbreg.b = pypris.observation.ravel()
    linbreg.PyPRIS_iter = "pris" + str(PRIS_iter)
    linbreg.stopping_loghistpercdelres_thres -= PRIS_iter
    linbreg.alpha = linbreg_ori.alpha / PRIS_iter
    print(" ")
    print(" ---------------- PRIS refinement # " + str(PRIS_iter) + " ------------------ ")
    print(" ")
    # recover
    linbreg.get_ready()
    linbreg.go()
