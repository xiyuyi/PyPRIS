import matplotlib.pyplot as plt
from PyPRIS import *

path = "G:/DH_localization/bgSCF1_mu5.0e+09/saved_objects" # specifie datafile position
PyPRIS_name = "PyPRIS_bgSCF1_mu5.0e+09_pris0_6678" # specify datafile name
PyPRIS_SensMx_name = "PyPRIS_bgSCF1_mu5.0e+09_pris0_SensingMx" # specify sensing matrix file name.

linbreg = loadCSSolver(path, PyPRIS_name, PyPRIS_SensMx_name)
linbreg.path_d = path
linbreg.debug = True
linbreg.debug_output(linbreg.it_count,'visualize')



v = linbreg.candidate_vis()
vis = v[:,:,:]
prj_ax0 = copy.deepcopy(np.mean(vis, axis=0))
prj_ax1 = copy.deepcopy(np.mean(vis, axis=1))
prj_ax2 = copy.deepcopy(np.mean(vis, axis=2).T)

patch = np.zeros((vis.shape[0],vis.shape[0]))

plt.figure(figsize=(3,3))
cat1 = np.concatenate([prj_ax0, prj_ax2], axis = 1)
cat2 = np.concatenate([prj_ax1, patch], axis = 1)
cat = np.concatenate([cat1, cat2], axis = 0)
plt.imshow(cat)
plt.savefig(
    '{}/PyPRIS__{}_{}_{}_plots_it{}.png'.format( linbreg.path_d, 'Proj_vies', linbreg.PyPRIS_name, linbreg.PyPRIS_iter, linbreg.it_count),
    dpi=300, figsize=(100, 80))

