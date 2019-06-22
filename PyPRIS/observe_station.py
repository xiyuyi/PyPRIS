from cmath import nan

import numpy as np
from scipy import interpolate

class ObserveStation:
    
    def __init__(self):
        self.biplane_observer = None
        self.monoplane_observer = None
        # create position for 4 different channel observers.
        self.channel_observer_1 = None
        self.channel_observer_2 = None
        self.channel_observer_3 = None
        self.channel_observer_4 = None
        self.observer_with_shift = None
        self.observer_with_CL = None
        self.dist = None

    class SingleObs:
        
        def __init__(self):
            self.psf = np.arange(1, 31).reshape(2, 3, 5)
            self.location = [0, 1, 2]
            self.imsize = [4, 4]
            self.psfz0 = 3
            self.debug = True
            self.edge_padding = False
            self.subpixel_shift = True

        def single_obs(self):
            """ 
            opts.psf   numpy.ndarray, 3D psf matrix, x-y should be odd numbers.  
            opts.location list, location coordinte for z-x-y
            opts.imsize list, x-y dimension of the observation image.
            opts.psfz0 list, with a single integer. 
            """

            if self.debug is True:
                print('---------------------- debug message ----------------------')
                print('Warning: Avoid rectangular PSF matrix in xy plane!')
                print('')

            if self.subpixel_shift is True:
                loc1 = np.int(np.floor(self.location[1]))
                loc2 = np.int(np.floor(self.location[2]))
                loc1_sps = self.location[1] - loc1  # the subpixel shift component
                loc2_sps = self.location[2] - loc2  # the subpixel shift component

            locz = self.location[0]  # location coordinates of depth.

            obs1w = self.imsize[0];  # first dimension size
            obs2w = self.imsize[1];  # second dimension size

            psfdim = self.psf.shape
            psf1hw = (psfdim[1] - 1) // 2  # half width of the psf matrix in dimension 1, excluding the center point.
            psf2hw = (psfdim[2] - 1) // 2  # half width of the psf matrix in dimension 2, excluding the center point.

            self.layerN = np.int(locz + self.psfz0)

            loc1sta = int(max(loc1 - psf1hw, 0))
            loc1end = int(min(loc1 + psf1hw + 1, obs1w))
            loc2sta = int(max(loc2 - psf2hw, 0))
            loc2end = int(min(loc2 + psf2hw + 1, obs2w))

            psf1sta = int(psf1hw - (loc1 - loc1sta));  # need to think about whether I should add 1 or not.
            psf1end = int(psf1hw + (loc1end - loc1));
            psf2sta = int(psf2hw - (loc2 - loc2sta));  # need to think about whether I should add 1 or not.
            psf2end = int(psf2hw + (loc2end - loc2));

            if self.debug is True:
                print('original loc1 is ' + str(self.location[1]))
                print('original loc2 is ' + str(self.location[2]))
                print('self.subpixel_shift is ' + str(self.subpixel_shift))
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
                locs1new = locs1[1:locs1.size - 1] - loc1_sps
                locs2new = locs2[1:locs2.size - 1] - loc2_sps
                # find the shifted stamp and assin        
                if self.debug is True:
                    print('-------------------------- Now interpoltion ')
                    print('     stamp shape is: ' + str(stamp.shape))
                    print('     locs1 shape is: ' + str(locs1.shape))
                    print('     locs2 shape is: ' + str(locs2.shape))
                    print('  locs1new shape is: ' + str(locs1new.shape))
                    print('  locs2new shape is: ' + str(locs2new.shape))

                # define interpolation function
                f = interpolate.interp2d(locs2, locs1, stamp, kind='cubic')
                stamp_sps = f(locs2new, locs1new)
                stamp[1:stamp.shape[0] - 1, 1:stamp.shape[1] - 1] = stamp_sps
                self.stamp = np.copy(stamp[psf1sta:psf1end, psf2sta:psf2end])
            else:
                self.stamp = np.copy(self.psf[self.layerN, psf1sta:psf1end, psf2sta:psf2end])

            if self.edge_padding is True:
                ev1 = self.psf[self.layerN, [1, psf1hw * 2], :].ravel()
                ev2 = self.psf[self.layerN, :, [1, psf2hw * 2]].ravel()
                edge_value = np.concatenate((ev1, ev2), axis=0).ravel().mean()

                self.obs = np.ones(self.imsize) * edge_value
                self.obs[loc1sta:loc1end, loc2sta:loc2end] = self.stamp
            else:
                self.obs = np.zeros(self.imsize)
                self.obs[loc1sta:loc1end, loc2sta:loc2end] = self.stamp

            if self.debug is True:
                print('shape of observation: ' + str(self.imsize))
                print('        shape of psf: ' + str(self.psf.shape))
                print('     z0 index in psf: ' + str(self.psfz0))
                print('              psf1hw: ' + str(psf1hw))
                print('              psf2hw: ' + str(psf2hw))
                print('')
                print('')

                self.psfwbox = np.copy(self.psf[self.layerN, :, :])
                c = max(self.psfwbox.ravel())
                self.psfwbox[psf1sta, psf2sta: psf2end] = c
                self.psfwbox[psf1end - 1, psf2sta: psf2end] = c
                self.psfwbox[psf1sta: psf1end, psf2sta] = c
                self.psfwbox[psf1sta: psf1end, psf2end - 1] = c

                self.obswbox = np.copy(self.obs)
                c = max(self.obswbox.ravel())
                self.obswbox[loc1sta, loc2sta: loc2end] = c
                self.obswbox[loc1end - 1, loc2sta: loc2end] = c
                self.obswbox[loc1sta: loc1end, loc2sta] = c
                self.obswbox[loc1sta: loc1end, loc2end - 1] = c

    def observe_biplane_prep(self, psf, single_image_size, deltaz_plane1, deltaz_plane2, psfz0,
                             observer_debugger, observer_edge_padding):
        # prepare an observer for biplane observation
        # this method will only be executed once in the preparation before calculating the sensing matrix.
        # this prep method provides an input window in the main pris script.
        self.biplane_observer = self.SingleObs() # this is the child class.
        self.biplane_observer.psf = np.copy(psf)
        self.biplane_observer.psfz0 = psfz0
        self.biplane_observer.debug = observer_debugger
        self.biplane_observer.imsize = single_image_size  # this should be the image size of the single plane observation
        self.biplane_observer.edge_padding = observer_edge_padding  # yes we want edge padding.
        self.biplane_observer.deltaz_1 = deltaz_plane1  # the observer now knows the plane shift for the first plane.
        self.biplane_observer.deltaz_2 = deltaz_plane2  # the observer now knows the plane shift for the second plane.

    def observe_biplane(self, loc):
        # take the biplane observation
        # this method will be passed into the sensing matrix generator, and
        # be executed iterative throughout the course of sensing matrix generation.
        #
        # get the depth positions of the two observation planes for the biplane_observer.
        loc1 = np.copy(loc);
        loc1[0] = np.copy(loc[0]) + self.biplane_observer.deltaz_1
        loc2 = np.copy(loc);
        loc2[0] = np.copy(loc[0]) + self.biplane_observer.deltaz_2

        # the biplane_observer now observe the first plane.
        self.biplane_observer.location = loc1  # focus at the position
        self.biplane_observer.single_obs()  # take the observation
        self.biplane_observer.observation1 = self.biplane_observer.obs.ravel()  # record this first observation

        # the biplane_observer now observe the second plane.
        self.biplane_observer.location = loc2  # focus at the position
        self.biplane_observer.single_obs()  # take the observation
        self.biplane_observer.observation2 = self.biplane_observer.obs.ravel()  # record this second observation

        # the biplane_observer now returns the observation
        observation = np.concatenate([self.biplane_observer.observation1, self.biplane_observer.observation2]).ravel()
        return observation

    def observe_monoplane_prep(self, psf, single_image_size, psfz0,
                             observer_debugger, observer_edge_padding):
        # prepare an observer for one-plane observation
        # this method will only be executed once in the preparation before calculating the sensing matrix.
        # this prep method provides an input window in the main pris script.
        self.monoplane_observer = self.SingleObs() # this is the child class.
        self.monoplane_observer.psf = np.copy(psf)
        self.monoplane_observer.psfz0 = psfz0
        self.monoplane_observer.debug = observer_debugger
        self.monoplane_observer.imsize = single_image_size  # this should be the image size of the single plane observation
        self.monoplane_observer.edge_padding = observer_edge_padding  # yes we want edge padding.

    def observe_monoplane(self, loc):
        # take the single plane observation
        # this method will be passed into the sensing matrix generator, and
        # be executed iterative throughout the course of sensing matrix generation.
        self.monoplane_observer.location = loc  # focus at the position
        self.monoplane_observer.single_obs()  # take the observation
        self.monoplane_observer.observation = self.monoplane_observer.obs.ravel()  # record this single plane observation
        return self.monoplane_observer.observation.ravel()

    def observe_two_channel_prep(self, psf1, size1, psf2, size2, psfz0,
                                 observer_debugger, observer_edge_padding):
        # prepare an observer for two channel observations
        # this method will only be executed once in the preparation before calculating the sensing matrix.
        # this prep method provides an input window in the main pris script.
        # we need to assign two different channel observers.
        # prepare the first channel observer
        self.channel_observer_1 = self.SingleObs()  # this is the child class.
        self.channel_observer_1.psf = np.copy(psf1)
        self.channel_observer_1.imsize = size1  # this should be the image size of the single plane observation
        self.channel_observer_1.psfz0 = psfz0
        self.channel_observer_1.debug = observer_debugger
        self.channel_observer_1.edge_padding = observer_edge_padding  # yes we want edge padding.

        # prepare the second channel observer
        self.channel_observer_2 = self.SingleObs()  # this is the child class.
        self.channel_observer_2.psf = np.copy(psf2)
        self.channel_observer_2.imsize = size2  # this should be the image size of the single plane observation
        self.channel_observer_2.psfz0 = psfz0
        self.channel_observer_2.debug = observer_debugger
        self.channel_observer_2.edge_padding = observer_edge_padding  # yes we want edge padding.

    def observe_two_channel(self, loc):
        # take the simultaneous two channel observation
        # this method will be passed into the sensing matrix generator, and
        # be executed iterative throughout the course of sensing matrix generation.
        # observe the first channel. channel_observer_1 is configured to the first channel
        self.channel_observer_1.location = loc  # focus at the position
        self.channel_observer_1.single_obs()  # take the observation
        self.channel_observer_1.observation = self.channel_observer_1.obs.ravel()  # record this first observation

        # observe the second channel (channel_observer_2 is configured to the second channel)
        self.channel_observer_2.location = loc  # focus at the position
        self.channel_observer_2.single_obs()  # take the observation
        self.channel_observer_2.observation = self.channel_observer_2.obs.ravel()  # record this first observation

        # Now returns the combined observations from two channels
        observation = np.concatenate([self.channel_observer_1.observation, self.channel_observer_2.observation]).ravel()
        return observation

    def observe_with_CL_prep(self, psf, single_iamge_size, psfz0,
                        dist_1_amplitd, dist_1_shift,
                        dist_2_amplitd, dist_2_shift,
                        observer_debugger, observer_edge_padding):
        self.observer_with_CL = self.SingleObs()  # this is the child class.
        self.observer_with_CL.psf = np.copy(psf)
        self.observer_with_CL.imsize = single_iamge_size  # this should be the image size of the single plane observation
        self.observer_with_CL.psfz0 = psfz0
        self.observer_with_CL.debug = observer_debugger
        self.observer_with_CL.edge_padding = observer_edge_padding
        self.dist.CL_A1 = dist_1_amplitd
        self.dist.CL_S1 = dist_1_shift
        self.dist.CL_A2 = dist_2_amplitd
        self.dist.CL_S2 = dist_2_shift
        return nan

    def observe_with_CL(self, loc):
        # update loc to incorporate field distortion and alignment
        # get an observer to observe with updated location coordiantes
        # self.channel_observer_2 = self.SingleObs()
        loc_shifted = copy.deepcopy(loc)
        loc_shifted[1] = loc[1]*self.dist.CL_A1 + self.dist.CL_S1# update location based on field distortion parameters.
        loc_shifted[2] = loc[2]*self.dist.CL_A2 + self.dist.CL_S2# update location based on field distortion parameters.
        self.observer_with_CL.location = loc_shifted  # focus at the position
        self.observer_with_CL.single_obs()  # take the observation
        self.observer_with_CL.observation = self.observer_with_CL.obs.ravel()  # record this first observation
        return self.observer_with_CL.observation

    def observe_with_shift_prep(self, psf, single_image_size, psfz0,
                                      shift_1, shift_2,
                                      observer_debugger, observer_edge_padding):
        # the grating doesn't cause field distortion, but causes translation based on the alignment
        # and the cropping of the image from the raw data.
        # such translation movement needs to be characterized from the experimental data and fed into this observer.
        self.observer_with_shift = self.SingleObs()  # this is the child class.
        self.observer_with_shift.psf = np.copy(psf)
        self.observer_with_shift.imsize = single_image_size  # this should be the image size of the single plane observation
        self.observer_with_shift.psfz0 = psfz0
        self.observer_with_shift.debug = observer_debugger
        self.observer_with_shift.edge_padding = observer_edge_padding
        self.dist.shift_1 = shift_1
        self.dist.shift_2 = shift_2
        self.observer_with_shift.debug = observer_debugger

    def observe_with_shift(self, loc):
        loc_shifted = copy.deepcopy(loc)
        loc_shifted[1] = loc[1] + self.dist.shift_1 # update location based on field translation parameters.
        loc_shifted[2] = loc[2] + self.dist.shift_2 # update location based on field translation parameters.
        self.observer_with_shift.location = loc_shifted  # focus at the position
        self.observer_with_shift.single_obs()  # take the observation
        self.observer_with_shift.observation = self.observer_with_shift.obs.ravel()  # record this first observation
        return self.observer_with_shift.observation

    def observe_with_CL_and_grating_prep(self,  psf_CL, imsize_CL, psfz0_CL,
                                                dist_1_amplitd, dist_1_shift,
                                                dist_2_amplitd, dist_2_shift,
                                                psf_shift, imsize_shift, psfz0_shift,
                                                shift_1, shift_2,
                                                observer_debugger, observer_edge_padding):

        # wrap the observation with the combination of CL and grating channels.
        self.observe_with_CL_prep(psf_CL, imsize_CL, psfz0_CL,
                        dist_1_amplitd, dist_1_shift,
                        dist_2_amplitd, dist_2_shift,
                        observer_debugger, observer_edge_padding)

        self.observe_with_shift_prep(psf_shift, imsize_shift, psfz0_shift,
                                      shift_1, shift_2,
                                      observer_debugger, observer_edge_padding)

    def observe_with_CL_and_grating(self, loc):
        obs_with_CL = self.observe_with_CL(loc)
        obs_with_shift = self.observe_with_shift(loc)
        observation = np.concatenate([obs_with_CL, obs_with_shift]).ravel()
        return observation