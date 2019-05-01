import matplotlib.pyplot as plt
from PyPRIS import *
f = list()
f.append("bgSCF0.8_mu1.0e+091")
f.append("bgSCF0.8_mu1.0e+101")
f.append("bgSCF0.8_mu1.0e+111")
f.append("bgSCF0.8_mu1.0e+121")
f.append("bgSCF0.8_mu8.0e-011")
f.append("bgSCF1.5_mu1.0e+091")
f.append("bgSCF1.5_mu1.0e+101")
f.append("bgSCF1.5_mu1.0e+111")
f.append("bgSCF1.5_mu1.0e+121")
f.append("bgSCF1.5_mu1.5e+001")
f.append("bgSCF1_mu1.0e+001")
f.append("bgSCF1_mu1.0e+091")
f.append("bgSCF1_mu1.0e+101")
f.append("bgSCF1_mu1.0e+111")
f.append("bgSCF1_mu1.0e+121")
f.append("bgSCF2_mu1.0e+091")
f.append("bgSCF2_mu1.0e+101")
f.append("bgSCF2_mu1.0e+111")
f.append("bgSCF2_mu1.0e+121")
f.append("bgSCF2_mu2.0e+001")



PyPRIS_SensMx_name = "PyPRIS_Demo_pris0_SensingMx" # specify sensing matrix file name.
itN = 0
for fitem in f:
    path = "G:/DH_localization/PyPRIS_tickets2/"+fitem+"/saved_objects"  # specifie datafile position
    for itN in np.arange(0,2200,200):
        PyPRIS_name = "PyPRIS_Demo_pris0_"+str(1+itN)  # specify datafile name
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
