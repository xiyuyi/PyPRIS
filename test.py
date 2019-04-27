from PyPRIS import *
print("hello world!")
import warnings
warnings.filterwarnings("ignore")

dpath = './PyPRIS_notebooks/test_data/psf.tif';
psf = io.imread(dpath)


# load observation
data_path = './PyPRIS_notebooks/test_data';
f = 'obsStack.tif';
fov = io.imread(data_path + '/' + f)
fov.shape
plt.figure(figsize=(15,5));
plt.subplot(131)
plt.imshow(fov[[1,3,5,6,10,15],:,:].reshape(6*fov.shape[1],fov.shape[2]).T)

# choose the blur
blur1 = fov[1,:,:];
blur2 = fov[6,:,:];


"""

100 nm step interval so there is a 500 nm separation. 20 layer separation between psf steps.

"""


blur = np.concatenate([blur1, blur2]).reshape(2*71,71);


# prepare a observer for pypris. [from one candidate, to one observation]
#  #  summon an observer from the ObserveStation (ObserveStation class) (Cheers~)
observer = ObserveStation()

#  #  configure this observer with biplane observation specifics
observer.observe_biplane_prep(psf, single_image_size = blur1.shape, \
                           deltaz_plane1 = -10, \
                           deltaz_plane2 = 10, \
                           psfz0 = 106, \
                           observer_debugger = False, \
                           observer_edge_padding = True)
# this observer will perform observation for the pypris object as seen below in the next section titled 'Initialize PyPRIS'.



# get a PyPRIS object:
pypris = PyPRIS()
pypris.observation = blur.ravel()
pypris.current_candidates_intervals = list([1, 4, 4])  # initialize the first intervals of neighboring candidate voxels.

# decide candidates for the first PRIS iteration. PRIS #0:
range_ind0 = np.arange(-36, 36, pypris.current_candidates_intervals[0]) # seems like for DH psf, use finnest intervals for the thickness dimension is better.
range_ind1 = np.arange(5,   66, pypris.current_candidates_intervals[1])
range_ind2 = np.arange(5,   66, pypris.current_candidates_intervals[2])

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
linbreg_ori = LinBreg(None)
linbreg_ori.debug = True
linbreg_ori.deep_debug = False
linbreg_ori.mu = 1000000000
linbreg_ori.obs_dim0 = blur.shape[0]
linbreg_ori.obs_dim1 = blur.shape[1]
linbreg_ori.alpha = 1e-8
linbreg_ori.maxit = 200
linbreg_ori.it_check_rem = 1
linbreg_ori.debug_it_int = 100
linbreg_ori.kick.ints = 10
linbreg_ori.kick.flag = True
linbreg_ori.kick.thres = 1e-3
linbreg_ori.save_obj_int = 100
linbreg_ori.save = True



# construct sensing matrix,
pypris.generate_sensing_mx()

# prepare the inner sparse recovery
pypris.current_A[:,len(pypris.current_candidates)] = 550
Anorm = pypris.current_A/500
linbreg = copy.deepcopy(linbreg_ori)
linbreg.candidate_coords = pypris.current_candidates
linbreg.candidate_intervals = pypris.current_candidates_intervals
linbreg.A = Anorm
linbreg.b = pypris.observation.ravel()
linbreg.id = 0

# recover
linbreg.get_ready()
linbreg.go()
