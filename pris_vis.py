import matplotlib.pyplot as plt
from PyPRIS import *
f = list()
f.append("bgSCF0.5_mu1.0e+10_alpha1.0e-09")
f.append("bgSCF0.5_mu1.0e+10_alpha1.0e-10")
f.append("bgSCF0.5_mu1.0e+10_alpha1.0e-11")
f.append("bgSCF0.8_mu1.0e+10_alpha1.0e-09")
f.append("bgSCF0.8_mu1.0e+10_alpha1.0e-10")
f.append("bgSCF0.8_mu1.0e+10_alpha1.0e-11")
f.append("bgSCF1.2_mu1.0e+10_alpha1.0e-09")
f.append("bgSCF1.2_mu1.0e+10_alpha1.0e-10")
f.append("bgSCF1.2_mu1.0e+10_alpha1.0e-11")
f.append("bgSCF1_mu1.0e+10_alpha1.0e-09")
f.append("bgSCF1_mu1.0e+10_alpha1.0e-10")
f.append("bgSCF1_mu1.0e+10_alpha1.0e-11")

#for ll in [1]:
for fitem in f:
    #fitem = f[-1]
    PyPRIS_SensMx_name = "PyPRIS_" + fitem + "_pris0_SensingMx"  # specify sensing matrix file name.
    path = "G:/DH_localization/PyPRIS_tickets_set2/"+fitem+"/saved_objects"  # specifie datafile position
    for itN in np.arange(0,62001,2000):
        for prisIter in np.arange(0,6):
            PyPRIS_name = "PyPRIS_" + fitem + "_pris"+str(prisIter)+"_" + str(1 + itN)  # specify datafile name
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
