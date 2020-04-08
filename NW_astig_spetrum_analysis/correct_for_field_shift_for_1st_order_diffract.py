# import tools
from nd2reader import ND2Reader
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from scipy.optimize import curve_fit
import copy
obj_fpath = 'G:\\DH_localization\\PyPRIS_CL_set1\\bgSCF12_mu1.0e+08_alpha1.0e-05_thres-17zrange-19to19_chosen\\saved_objects'
linbreg_fpath = obj_fpath+'\\PyPRIS_bgSCF12_mu1.0e+08_alpha1.0e-05_thres-17zrange-19to19_pris6_80000.file'
pypris_fpath = obj_fpath+'\\PyPRIS_pris6.file'

'''

load the linbreg object from the 0th order result

'''

with open(linbreg_fpath,'rb') as f:
    linb0 = pickle.load(f)

with open(pypris_fpath,'rb') as f:
    pypris0 = pickle.load(f)

'define a new observer for the pypris object'