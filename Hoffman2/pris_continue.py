#!/u/home/x/xiyuyi/.conda/envs/PyPRIS_env/bin/python
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
print("PyPRIS Go!")
import warnings
warnings.filterwarnings("ignore")

# take a PRIS recovery order
ticket = get_ticket('./Go.pris_ticket')

# prep the order according to the ticket
psf = io.imread(ticket.psf_path)
blur1 = io.imread(ticket.plane1_path)
blur2 = io.imread(ticket.plane2_path)
blur = np.concatenate([blur1, blur2],axis=0)

# prepare a observer for pypris. [from one candidate, to one observation]
# summon an observer from the ObserveStation (ObserveStation class) (Cheers~)
observer = ObserveStation()
"configure this observer with biplane observation specifics"
observer.observe_biplane_prep(psf, single_image_size = blur1.shape,
                           deltaz_plane1 = ticket.plane1_dz,
                           deltaz_plane2 = ticket.plane2_dz,
                           psfz0 = ticket.psfz0,
                           observer_debugger = ticket.observer_debugger,
                           observer_edge_padding = ticket.observer_edge_padding)


# this observer will perform observation for the pypris object as seen below in the next section titled 'Initialize PyPRIS'.
# get a PyPRIS object:
pypris = PyPRIS()
pypris.observation = blur.ravel()
pypris.current_candidates_intervals = ticket.init_candidates_intervals  # initialize the first intervals of neighboring candidate voxels.

# decide candidates for the first PRIS iteration. PRIS #0:
range_ind0 = np.arange(ticket.init_ax0_range[0], ticket.init_ax0_range[1], pypris.current_candidates_intervals[0]) # seems like for DH psf, use finnest intervals for the thickness dimension is better.
range_ind1 = np.arange(ticket.init_ax1_range[0], ticket.init_ax1_range[1], pypris.current_candidates_intervals[1])
range_ind2 = np.arange(ticket.init_ax2_range[0], ticket.init_ax2_range[1], pypris.current_candidates_intervals[2])

# initialize the first pool of candidates with the current candidate intervals.
pypris.current_candidates = list()  # current pool of candidates.
for i0 in range_ind0:
    for i1 in range_ind1:
        for i2 in range_ind2:
            pypris.current_candidates.append([i0,i1,i2])

# set a check mark for pypris:
pypris.set_check_mark()

# put the configured observer to the pypris to perform the tasks of 'observe_biplane'
# when a pypris object observes, it is through the observer and it is
# done with the observer's observe_biplane skill (skill = method).
pypris.observe = observer.observe_biplane
#pypris.show_attributes()



# condition a main linbreg object that carries the basic configurations.
# for each PRIS iteration, a new linbreg object will be cloned from this one for its own task.
# the purpose here is just to simplify the parameter configuration process.
linbreg_ori = ticket.linbreg_alpha
linbreg_ori.obs_dim0 = blur.shape[0]
linbreg_ori.obs_dim1 = blur.shape[1]
linbreg_ori.PyPRIS_name = ticket.name

# construct sensing matrix,
pypris.generate_sensing_mx()
# prepare the inner sparse recovery
pypris.current_A = pypris.current_A/ticket.psf_norm_factor
c = np.dot(pypris.observation.ravel(),pypris.current_A)
c1 = np.min(c[0:-2])
c2 = np.max(c[0:-2])
bgv = c1 + (c2-c1)*ticket.bg_scaling_coef
pypris.current_A[:,-1] = bgv/np.sum(pypris.observation.ravel()) # the last column in A corresponds to the background component.

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
print("PRIS continue! now load a linbreg object from last round")
path = "/u/scratch/x/xiyuyi/PyPRIS_tickets_set4/bgSCF1.5_mu1.0e+10_alpha1.0e-09/saved_objects"  # specifie datafile position
PyPRIS_name = "PyPRIS_bgSCF1.5_mu1.0e+10_alpha1.0e-09_pris0_80001"  # specify datafile name
PyPRIS_SensMx_name = "PyPRIS_bgSCF1.5_mu1.0e+10_alpha1.0e-09_pris0_SensingMx"  # specify datafile name
linbreg = loadCSSolver(path, PyPRIS_name, PyPRIS_SensMx_name)

# construct sensing matrix,
for PRIS_iter in np.arange(1,ticket.PRIS_iter_end):
    pypris.prep_for_new_refinement()
    pypris.refine_candidates(linbreg)
    pypris.generate_sensing_mx()
    # prepare the inner sparse recovery
    pypris.current_A = pypris.current_A / ticket.psf_norm_factor
    c = np.dot(pypris.observation.ravel(), pypris.current_A)
    c1 = np.min(c[0:-2])
    c2 = np.max(c[0:-2])
    bgv = c1 + (c2 - c1) * ticket.bg_scaling_coef
    pypris.current_A[:, -1] = bgv / np.sum(pypris.observation.ravel())  # the last column in A corresponds to the background component.
    linbreg = copy.deepcopy(linbreg_ori)
    linbreg.candidate_coords = pypris.current_candidates
    linbreg.candidate_intervals = pypris.current_candidates_intervals
    linbreg.A = pypris.current_A
    linbreg.b = pypris.observation.ravel()
    linbreg.PyPRIS_iter = "pris"+ str(PRIS_iter)
    linbreg.stopping_loghistpercdelres_thres -= PRIS_iter*2
    print("---------------- PRIS refinement #" + str(PRIS_iter) + " ------------------")
    # recover
    linbreg.get_ready()
    linbreg.go()
