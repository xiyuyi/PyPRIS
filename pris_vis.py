import os
import matplotlib.pyplot as plt
from PyPRIS import *
import copy

"""

create png files for all the results in the specified directory and sub directories

"""
# put in the file folder where you want all the linbreg objects to be visualized.
path = 'G:\\DH_localization\\PyPRIS_tickets_F2_set1'


files = []
filenames = []
filepaths = []
# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for file in f:
        if file.endswith('1.file'):
            files.append(os.path.join(r, file))
            filenames.append(file[0:-5])
            filepaths.append(r)

for f, fname, fpath in zip(files, filenames, filepaths):
    print(f)
    print(fname)
    print(fpath)
    print("--------------")

for path, PyPRIS_name in zip(filepaths, filenames):
    #PyPRIS_name =  "PyPRIS_" + fitem + "_pris"+str(prisIter)+"_" + str(1 + itN)  # specify datafile name
    tp = PyPRIS_name.split('_')
    PyPRIS_SensMx_name = "_".join(tp[0:-1])+"_SensingMx"  # specify datafile name
    print(PyPRIS_name)
    try:
        linbreg = loadCSSolver(path, PyPRIS_name, PyPRIS_SensMx_name)
        linbreg.path_d = path
        linbreg.debug = True
        try:
            linbreg.debug_output(linbreg.it_count,'visualize')
        except:
            pass
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
    except:
        pass
