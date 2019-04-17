import numpy as np
import copy
from scipy import interpolate

class SingleObs:
    def __init__(self):
        self.psf = np.arange(1,31).reshape(2,3,5)
        self.location = [0,1,2]
        self.imsize=[4,4]
        self.psfz0 = 3
        self.debug = False
        self.edge_padding = False
        self.subpixel_shift = False
        
    def single_obs(self):
        """ 
        opts.psf   numpy.ndarray, 3D psf matrix, x-y should be odd numbers.  
        opts.location list, location coordinte for z-x-y
        opts.imsize list, x-y dimension of the observation image.
        opts.psfz0 list, with a single integer. 
        """

        if self.debug is True:
            print('----------------- debug message -----------------')
            print('Warning: Avoid rectangular PSF matrix in xy plane!')
            print('')
            
        if self.subpixel_shift is True:
            loc1 = np.int(np.floor(self.location[1]))
            loc2 = np.int(np.floor(self.location[2]))
            loc1_sps = self.location[1] - loc1 # the subpixel shift component
            loc2_sps = self.location[2] - loc2 # the subpixel shift component
            if loc1_sps == 0 and loc2_sps == 0:
                self.subpixel_shift = False
                
        else:
            loc1 = self.location[1] # first dimension location coordiantes. [inner]
            loc2 = self.location[2] # second dimension location coordiantes. [further inner]

        locz = self.location[0] # location coordinates of depth.
        
        obs1w = self.imsize[0]; # first dimension size
        obs2w = self.imsize[1]; # second dimension size

        psfdim = self.psf.shape
        psf1hw = (psfdim[1] - 1) // 2 # half width of the psf matrix in dimension 1, excluding the center point.
        psf2hw = (psfdim[2] - 1) // 2 # half width of the psf matrix in dimension 2, excluding the center point.

        self.layerN = np.int(locz + self.psfz0)

        loc1sta = int(max(loc1 - psf1hw, 0))
        loc1end = int(min(loc1 + psf1hw + 1, obs1w))
        loc2sta = int(max(loc2 - psf2hw, 0))
        loc2end = int(min(loc2 + psf2hw + 1, obs2w))

        psf1sta = int(psf1hw - (loc1 - loc1sta)); # need to think about whether I should add 1 or not.
        psf1end = int(psf1hw + (loc1end - loc1)); 
        psf2sta = int(psf2hw - (loc2 - loc2sta)); # need to think about whether I should add 1 or not.
        psf2end = int(psf2hw + (loc2end - loc2)); 
        
        if self.debug is True:
            print('original loc1 is ' + str(self.location[1]))
            print('original loc2 is ' + str(self.location[2]))
            print('loc2 is ' + str(loc2))
            print('loc1 is ' + str(loc1))
            print('loc2 is ' + str(loc2))
            print('loc1 excess is ' + str(loc1_sps))
            print('loc2 excess is ' + str(loc2_sps))
            print(' location coordiante: locz=' + str(locz) + ', loc1=' + str(loc1) + ', loc2=' + str(loc2))
            print('')
            print('layerN = ' + str(self.layerN))
            print('loc1sta = ' + str(loc1sta) + '; loc1end = ' + str(loc1end))
            print('loc2sta = ' + str(loc2sta) + '; loc2end = ' + str(loc2end))
            print('psf1sta = ' + str(psf1sta) + '; psf1end = ' + str(psf1end))
            print('psf2sta = ' + str(psf2sta) + '; psf2end = ' + str(psf2end))
            
        if self.subpixel_shift is True:
            # get the stamp with only the integer part
            stamp = np.copy(self.psf[self.layerN, :, :])
            # up perform sub-pixel interpolation to include the sub-pixel shifts in the stamp.
                
            locs1 = np.arange(0, stamp.shape[0])
            locs2 = np.arange(0, stamp.shape[1])
            # find interpolation coordinates based on the sub-pixel shifts
            locs1new = locs1[1:locs1.size-1] - loc1_sps
            locs2new = locs2[1:locs2.size-1] - loc2_sps
            # find the shifted stamp and assin        
            if self.debug is True:
                print('--------------------------Now interpoltion ')
                print('     stamp shape is: ' + str(stamp.shape))
                print('     locs1 shape is: ' + str(locs1.shape))
                print('     locs2 shape is: ' + str(locs2.shape))
                print('  locs1new shape is: ' + str(locs1new.shape))
                print('  locs2new shape is: ' + str(locs2new.shape))
                
            
            # define interpolation function
            f = interpolate.interp2d(locs2, locs1, stamp, kind='cubic')
            stamp_sps = f(locs2new, locs1new)
            stamp[1:stamp.shape[0]-1,1:stamp.shape[1]-1] = stamp_sps;
            
            self.stamp = np.copy(stamp[psf1sta:psf1end, psf2sta:psf2end])
        else:
            self.stamp = np.copy(self.psf[self.layerN, psf1sta:psf1end, psf2sta:psf2end])
        
        if self.edge_padding is True:
            ev1 = self.psf[self.layerN,[1,psf1hw*2],:].ravel()
            ev2 = self.psf[self.layerN,:,[1,psf2hw*2]].ravel()
            edge_value = np.concatenate((ev1,ev2),axis=0).ravel().mean()

            self.obs = np.ones(self.imsize)*edge_value 
            self.obs[loc1sta:loc1end,loc2sta:loc2end] = self.stamp
        else:
            self.obs = np.zeros(self.imsize)
            self.obs[loc1sta:loc1end,loc2sta:loc2end] = self.stamp
        
        if self.debug is True:
            print('shape of observation: ' + str(self.imsize))
            print('        shape of psf: ' + str(self.psf.shape))
            print('     z0 index in psf: ' + str(self.psfz0))
            print('     z0 index in psf: ' + str(self.psfz0))
            print('              psf1hw: ' + str(psf1hw))
            print('              psf2hw: ' + str(psf2hw))
            print('')
            print('')
            
            self.psfwbox = np.copy(self.psf[self.layerN, :, :])
            c = max(self.psfwbox.ravel())
            self.psfwbox[psf1sta,     psf2sta : psf2end] = c
            self.psfwbox[psf1end - 1, psf2sta : psf2end] = c
            self.psfwbox[psf1sta : psf1end,     psf2sta] = c
            self.psfwbox[psf1sta : psf1end, psf2end - 1] = c
            
            self.obswbox = np.copy(self.obs)
            c = max(self.obswbox.ravel())
            self.obswbox[loc1sta,     loc2sta : loc2end] = c
            self.obswbox[loc1end - 1, loc2sta : loc2end] = c
            self.obswbox[loc1sta : loc1end,     loc2sta] = c
            self.obswbox[loc1sta : loc1end, loc2end - 1] = c
            
        
class PyPRIS:
    def __init__(self):
        self.name = 'PRIS object'
        self.A = np.ndarray(0) # sensing matrix
        self.current_ReF = 0 # refinement factor
        self.candidate_coordinates = list() # coordinates of candidates
        self.candidate_3D_space_Ranges = dict() # the range of the coordinates of candidates
        self.blur = np.ndarray() # observation
        self.positivity = True # positivity constraint

        
