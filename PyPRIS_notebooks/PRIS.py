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
        self.candidate_3D_space_Ranges = list() # the range of the coordinates of candidates
        self.positivity = True # positivity constraint
    def test(self):
        print('this is a test')
        
class LinBreg:
    # Authors: Xiyu Yi, Xingjia Wang @ UCLA, 2019.
    # PI: Shimon Weiss, Department of Chemistry and Biochemistry, UCLA.
    import time
    def __init__(self, A, x, b):
        # solve for x from Ax = b.
        self.A = A # sensing matrix.
        self.x = x # coefficient vector for the pool of candidates.
        self.b = b # observation vector.
        self.mu = np.mean(self.b.ravel())  # shrinkage threshold.
        self.er = np.zeros(self.b.shape)
        self.erpj = np.zeros(self.x.shape)
        self.cumerr = np.zeros(self.x.shape)
        self.recb = np.zeros(self.b.shape)
        self.flag_stop = False # flag to stop optimization iteration.
        self.maxit = 2 # maximum iteration steps.
        self.debug_it_int = 1
        self.flag_positivity = True
        self.it_check_rem = 1
        self.iterations = list()
        self.hist_res = list()
        self.hist_resDrop = list()
        
        self.bg = list()
        self.alpha = 1
        self.debug = False
        self.deep_debug = False
        self.kicking_flag = False
        self.kicking_ints = 10 # number of iterations between kicking evaluation.
        self.kicking_reference = self.x
        self.kicking_option = False
        self.kicking_thres = 1e-10
        self.kicking_refnorm = np.max(self.x)
        self.hist_kicking_refnorm = list()
        self.hist_kicking_eval_counts = list()
        self.stepsize = np.ones(self.x.shape) # stepsize.     

    def getready(self):
        if self.debug is True:
            import os
            # define the name of the directory to be created.
            path = "../../PyPRIS_Scratch/debug_output"
            try:  
                os.mkdir(path)
            except OSError:  
                print ("Creation of the directory %s failed" % path)
            else:  
                print ("Successfully created the directory %s " % path)
                

    def shrink(self,sk):
        sk[np.where((sk >= -self.mu) * (sk <= self.mu))] = 0
        sk[np.where(sk > self.mu)] -= self.mu
        sk[np.where(sk < -self.mu)] += self.mu
        return sk
    
    def kicking_evaluation(self):
        self.kicking_refnorm = np.linalg.norm(self.x-self.kicking_reference)
        if self.kicking_refnorm < self.kicking_thres:
            # flip the kicking flag to "True" with positive evaluation
            self.kicking_flag = True
        else:
            # flip the kicking flag to "False" with negative evaluation.
            self.kicking_flag = False
            
        self.hist_kicking_refnorm.append(self.kicking_refnorm)
        self.hist_kicking_eval_counts.append([it_count, self.kicking_flag])

    def kicking_go(self):
        # execute kicking
        # kicking execution only modifies the domains of the step size 
        # and split domains in a binary way.
        i0 = np.where(self.x == 0) # zero entries on x [x entries where there kicking is in demand]
        i1 = np.where(self.x != 0) # none zero entries in x. [x entries where there is a value, no need for kicking]
        si = (self.mu *np.sign(self.erpj[i0]) - self.x[i0]) / self.erpj[i0] # stepsie for entries that needs kicking
        self.stepsize[i0] = np.min(si)  
        self.stepsize[i1] = 1

        # reset kicking_flag to False and wait for the flip 
        # from the next positive kicking evaluation
        self.kicking_flag = False 

        # update kicking_reference for follow-up kicking evaluation
        self.kicking_reference = copy.deepcopy(self.x)

    def go(self):
        t1 = time.time()
        it_count = 0
        self.hist_res.append(0)
        self.hist_resDrop.append(0)
        self.iterations.append(it_count)
        self.bg.append(self.x[self.x.size-1]) 
        # main linearized bregman iteration with kicking option.
        while self.flag_stop is False:
            # incrementation of the iteration number.
            it_count += 1
            
            # calculate distance (error)
            self.recb = np.dot(self.A, self.x)
            self.er = self.b - self.recb 
            if self.deep_debug is True: self.debug_output(it_count, appstr = '_a_er_updated')
            
            # perform back projection of the error ('adding the errors back').
            self.erpj = np.dot(self.er, self.A)  
            if self.deep_debug is True: self.debug_output(it_count, appstr = '_b_erpj_updated')
            
            # check if kicking is needed 
            #
            # "Kicking" rescales the backprojected error (self.erpj) with two different stepsizes
            # we'll have stepsize > 1 for kicking area, and stepsize = 1 for non-kicking area. 
            # kicking boosts the tip of the cumulated backprojected errors towards the shrinkage 
            # threshold.
            # In this implementation, the effect of kicking is realized throught a modified. 
            # distribution of stepsizes (self.stepsize). 
            if np.remainder(it_count, self.kicking_ints) == 0:
                self.kicking_evaluation()
                # kick if we get a positive kicking ealuation.
                if self.kicking_flag is True: self.kicking_go()
                    
            # get the acumulation of the back projected error.
            self.cumerr += self.erpj * self.stepsize
            if self.deep_debug is True: self.debug_output(it_count, appstr = '_c_cumerr_updated')
            
            # perform positivity constraint:
            if self.flag_positivity is True: self.cumerr[np.where(self.cumerr < 0)] = 0
            if self.deep_debug is True: self.debug_output(it_count, appstr = '_d_positivity_updated')

            # shrinkage to update the candidate coefficients.
            self.x = copy.deepcopy(self.cumerr)
            if self.deep_debug is True: self.debug_output(it_count, appstr = '_e_x_copied')
            self.x = self.alpha * self.shrink(self.x)   
            if self.deep_debug is True: self.debug_output(it_count, appstr = '_f_x_updated')
            
            # decide on the termination of iterations.
            if it_count > self.maxit: self.flag_stop = True
            
            # update the quantities for status tracking purposes.
            self.track_status(it_count, self.er)
            
            # check intermediate outputs. (Valid under debug mode).
            self.debug_output(it_count, appstr = '_g_track_status_updated')

    # Generate intermediate output under debug mode.
    def debug_output(self, it_count, appstr):
        if self.debug is True:
            if np.remainder(it_count, self.debug_it_int) == self.it_check_rem:
                print('intermediate output it#'+ str(it_count))
                temp = np.mean(self.x[0:self.x.size-1].reshape(range_ind0.size,range_ind1.size,range_ind2.size),axis=0)
                nrow = 3
                ncol = 4
                plt.figure(figsize=(11,7))
                plt.subplot(nrow, ncol, 1)
                plt.imshow(temp)
                t = plt.title('XY-plane projection')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 2)
                plt.plot(self.iterations, np.log(self.hist_res), '.')
                t = plt.title('Log(L2(res))')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 3)
                plt.plot(self.er.ravel(), '.')
                t = plt.title('err')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 4)
                plt.text(0,0.8,'mu: '+str(np.floor(self.mu)),fontsize=16)
                plt.text(0,0.6,'stepsize: '+str(np.floor(self.stepsize)),fontsize=16)
                
                plt.axis('off')

                plt.subplot(nrow, ncol, 5)
                plt.plot(self.iterations, self.bg, '.')
                t = plt.title('Background')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 6)
                plt.plot(self.x.ravel())
                t = plt.title('coefficients (x)')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 7)
                plt.plot(self.cumerr.ravel(), '.')
                t = plt.title('cum-err')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 8)
                plt.plot(self.cumerr.ravel(), '.')
                plt.plot([0,len(self.cumerr.ravel())],[self.mu,self.mu],'r')
                t = plt.title('cum-err and mu')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 9)
                plt.plot(self.erpj.ravel(), '.')
                t = plt.title('erorr back projection (erpj)')
                t.set_position([.5, 1.15])
                
                plt.subplot(nrow, ncol, 10)
                plt.plot(self.b.ravel(), '.')
                t = plt.title('input obsrvation plot')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 11)
                plt.plot(self.recb.ravel(), '.')
                t = plt.title('recovered obsrvation plot')
                t.set_position([.5, 1.15])

                plt.tight_layout(rect=[0, 0.04, 1, 0.9])
                plt.subplots_adjust(top=0.85, left = 0.1)
                plt.savefig('../../PyPRIS_Scratch/debug_output/plots_it' + str(it_count) + appstr +'.png', dpi=300, figsize=(100,80))
                plt.close()

    def track_status(self, it_count, er):
        self.hist_res.append(np.linalg.norm(er))
        self.hist_resDrop.append((self.hist_res[it_count] - self.hist_res[it_count-1])/self.hist_res[it_count-1])
        self.iterations.append(it_count)
        self.bg.append(self.x[self.x.size-1]) 