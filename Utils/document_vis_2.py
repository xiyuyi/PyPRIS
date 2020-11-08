import imageio
# visualize the 4 field of views, as the blur inputs.
import matplotlib
matplotlib.use('agg')
import sys
sys.path.append("G:\\DH_localization\\PyPRIS")
from PyPRIS import *
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.misc
pixel_size = 62.75 # pixel size is 62.75 nanometers. notes: Measure pixel size #79   link: https://github.com/xiyuyi/PyPRIS/issues/79
pathAll= "G:\DH_localization\Experiments_Phase2_afterSPIE\April2019_BoulderTrip_analysis_choose_first_FOV_initialTest"

linbreg_paths = list()
linbreg_paths.append(pathAll + "PyPRIS_bgSCF1.5_mu1.0e+05_alpha1.0e-11_pris0_36838.file")
linbreg_paths.append(pathAll + "PyPRIS_bgSCF1.5_mu1.0e+05_alpha1.0e-11_pris1_75208.file")
linbreg_paths.append(pathAll + "PyPRIS_bgSCF1.5_mu1.0e+05_alpha1.0e-11_pris2_400000.file")
linbreg_paths.append(pathAll + "PyPRIS_bgSCF1.5_mu1.0e+05_alpha1.0e-11_pris3_400000.file")
output_path= ticket_path

prj_ax0s = []
prj_ax1s = []
prj_ax2s = []
blist = []
name_heads = []
candidate_intervals = []
for path in linbreg_paths:
    with open(path,'rb') as f:
        linbreg=pickle.load(f)
    name_heads.append(path[-17:-5])
    vis = linbreg.candidate_vis()
    prj_ax0s.append(copy.deepcopy(np.mean(vis, axis=0)))
    prj_ax1s.append(copy.deepcopy(np.mean(vis, axis=1)))
    prj_ax2s.append(copy.deepcopy(np.mean(vis, axis=2).T))
    blist.append(linbreg.b)
    candidate_intervals.append(linbreg.candidate_intervals)

for name_head, prj_ax0, prj_ax1, prj_ax2, b, interv in zip(name_heads, prj_ax0s, prj_ax1s, prj_ax2s, blist, candidate_intervals):
    this_path = output_path + '/proj_ax0_' + name_head + '_pixelsize'+str(pixel_size*interv[1]) + '.png'
    im = np.kron(prj_ax0, np.ones((100, 100)))
    #draw a scale bar of length 800 nm
    pcounts = np.round(300/(interv[2]*pixel_size/100))
    ysta = np.int(im.shape[0] * 0.92)
    yend = np.int(im.shape[0] * 0.96)
    xend = np.int(im.shape[1] * 0.95)
    im[ysta:yend,xend-np.int(pcounts):xend]=np.max(im.ravel())
    scipy.misc.toimage(im, cmin=0.0, cmax=np.max(prj_ax0.ravel())).save(this_path)
    #
    #
    this_path = output_path + '/proj_ax1' + name_head + '.png'
    im = np.kron(prj_ax1, np.ones((100, 100))).T
    # draw a scale bar of length 800 nm
    pcounts = np.round(300 / (interv[0] * z_step_size / 100))
    ysta = np.int(im.shape[0] * 0.92)
    yend = np.int(im.shape[0] * 0.96)
    xend = np.int(im.shape[1] * 0.95)
    im[ysta:yend,xend-np.int(pcounts):xend]=np.max(im.ravel())
    scipy.misc.toimage(im.T, cmin=0.0, cmax=np.max(prj_ax1.ravel())).save(this_path)
    #
    this_path = output_path + '/proj_ax2' + name_head + '.png'
    im = np.kron(prj_ax2.T, np.ones((100, 100))).T
    # draw a scale bar of length 800 nm
    pcounts =  np.round(300 /(interv[0] * z_step_size / 100))
    ysta = np.int(im.shape[0] * 0.92)
    yend = np.int(im.shape[0] * 0.96)
    xend = np.int(im.shape[1] * 0.95)
    im[ysta:yend, xend - np.int(pcounts):xend] = np.max(im.ravel())
    scipy.misc.toimage(im, cmin=0.0, cmax=np.max(prj_ax2.ravel())).save(this_path)
    #
    this_path = output_path + '/b' + name_head + '.png'
    d1 = 81
    d2 = 81

    im = np.kron(np.resize(b,(d1*2,d2)), np.ones((100, 100))).T
    # draw a scale bar of length 800 nm
    pcounts = np.round(300 / (1 * pixel_size / 100))
    ysta = np.int(im.shape[0] * 0.92)
    yend = np.int(im.shape[0] * 0.96)
    xend = np.int(im.shape[1] * 0.95)
    im[ysta:yend, xend - np.int(pcounts):xend] = np.max(im.ravel())
    im[10:11 + np.int(pcounts), 3] = np.max(im.ravel())
    scipy.misc.toimage(im.T, cmin=0.0, cmax=np.max(b.ravel())).save(this_path)



