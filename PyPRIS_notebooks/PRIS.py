import numpy as np
from scipy import interpolate
import time
import copy
import pickle
from matplotlib import pyplot as plt
import matplotlib

#  Authors: Xiyu Yi, Xingjia Wang @ UCLA, 2019.
#  PI: Shimon Weiss, Department of Chemistry and Biochemistry, UCLA.


class ObserveStation:
    def __init__(self):
        self.biplane_observer = None

    class SingleObs:
        def __init__(self):
            self.psf = np.arange(1, 31).reshape(2, 3, 5)
            self.location = [0, 1, 2]
            self.imsize = [4, 4]
            self.psfz0 = 3
            self.debug = False
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
                print('----------------- debug message -----------------')
                print('Warning: Avoid rectangular PSF matrix in xy plane!')
                print('')

            if self.subpixel_shift is True:
                loc1 = np.int(np.floor(self.location[1]))
                loc2 = np.int(np.floor(self.location[2]))
                loc1_sps = self.location[1] - loc1  # the subpixel shift component
                loc2_sps = self.location[2] - loc2  # the subpixel shift component
                if loc1_sps == 0 and loc2_sps == 0:
                    self.subpixel_shift = False

            else:
                loc1 = self.location[1]  # first dimension location coordiantes. [inner]
                loc2 = self.location[2]  # second dimension location coordiantes. [further inner]

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
                locs1new = locs1[1:locs1.size - 1] - loc1_sps
                locs2new = locs2[1:locs2.size - 1] - loc2_sps
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
                stamp[1:stamp.shape[0] - 1, 1:stamp.shape[1] - 1] = stamp_sps;

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

    def observe_biplane_prep(self, psf, single_image_size, deltaz_plane1=-10, deltaz_plane2=10, psfz0=106,\
                             observer_debugger=False, observer_edge_padding=True):
        # prepare an observer for biplane observation
        # this method will only be executed once in the preparation before calculating the sensing matrix.
        self.biplane_observer = self.SingleObs() # this is the parent class.
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
        blur = np.concatenate([self.biplane_observer.observation1, self.biplane_observer.observation2]).ravel()
        return blur


class PyPRIS:
    def __init__(self):
        self.name = 'PyPRIS object'
        self.positivity = True  # positivity constraint.
        self.observation = np.ndarray(0)  # this should be the observation vector

        self.current_relReF = list([1, 2, 2])    # current relative refinement factor. This is the relative refinement to be or have been performed
        #                               for this round of pris as compared to the last round of pris.
        self.current_PRIS_ItN = list([0])  # current PRIS iteration count, starting from 0.
        self.current_A = np.ndarray(0)  # current sensing matrix.
        self.current_candidates = list()  # current pool of candidates.
        self.current_candidates_intervals = list()  # current intervals of neighboring candidate voxels.
        self.current_check_mark_id = 0

        self.hist_candidates = list()  # keep a record of the full history.
        self.hist_candidates_intervals = list()  # keep a record of the full history.
        self.hist_PRIS_ItN = list()
        self.hist_check_mark_id = list()  # checkmark ID. ascending each time after you make a check_mark

        self.observator = None  # this should be a function that needs to be defined.

    def prep_for_new_refinement(self):
        self.hist_candidates.append(copy.deepcopy(self.current_candidates))
        self.hist_candidates_intervals.append(copy.deepcopy(self.current_candidates_intervals))
        self.hist_PRIS_ItN.append(copy.deepcopy(self.current_PRIS_ItN))
        self.set_check_mark()

    def refine_candidates(self, linbreg):
    # this will take the current candidates and the result in the linbreg object,
    # and generate refined pool of candidates.
        # create a check mark.
        self.set_check_mark()

        # Get the non_zero_coordinates from the existing linbreg results.
        non_zero_inds = np.argwhere(linbreg.x[0:len(linbreg.x) - 1] > 0)
        non_zero_coordinates = [self.current_candidates[i] for i in list(non_zero_inds.ravel())]

        self.current_candidates_intervals = copy.deepcopy(
            [pre / ref for pre, ref in zip(self.current_candidates_intervals, self.current_relReF)]
        )
        current_interval = self.current_candidates_intervals
    # get new coordinates with 2-fold refinement.
        new_coords = list()
        for i in non_zero_coordinates:
            extra_coords = [[i[0], i[1] - current_interval[1] / 2 * 3, i[2] - current_interval[2] / 2 * 3], \
                            [i[0], i[1] - current_interval[1] / 2 * 3, i[2] - current_interval[2] / 2], \
                            [i[0], i[1] - current_interval[1] / 2 * 3, i[2] + current_interval[2] / 2], \
                            [i[0], i[1] - current_interval[1] / 2 * 3, i[2] + current_interval[2] / 2 * 3], \
                            [i[0], i[1] - current_interval[1] / 2, i[2] - current_interval[2] / 2 * 3], \
                            [i[0], i[1] - current_interval[1] / 2, i[2] - current_interval[2] / 2], \
                            [i[0], i[1] - current_interval[1] / 2, i[2] + current_interval[2] / 2], \
                            [i[0], i[1] - current_interval[1] / 2, i[2] + current_interval[2] / 2 * 3], \
                            [i[0], i[1] + current_interval[1] / 2 * 3, i[2] - current_interval[2] / 2 * 3], \
                            [i[0], i[1] + current_interval[1] / 2 * 3, i[2] - current_interval[2] / 2], \
                            [i[0], i[1] + current_interval[1] / 2 * 3, i[2] + current_interval[2] / 2], \
                            [i[0], i[1] + current_interval[1] / 2 * 3, i[2] + current_interval[2] / 2 * 3], \
                            [i[0], i[1] + current_interval[1] / 2, i[2] - current_interval[2] / 2 * 3], \
                            [i[0], i[1] + current_interval[1] / 2, i[2] - current_interval[2] / 2], \
                            [i[0], i[1] + current_interval[1] / 2, i[2] + current_interval[2] / 2], \
                            [i[0], i[1] + current_interval[1] / 2, i[2] + current_interval[2] / 2 * 3]]
            for i1 in extra_coords:
                if i1 not in new_coords:
                    new_coords.append(i1)

        self.current_candidates = copy.deepcopy(new_coords)

        # set a check mark for tracking purposes.
        self.set_check_mark()

    def generate_sensing_mx(self):
        self.current_A = np.ndarray([len(self.observation), len(self.current_candidates) + 1])
        for count, loc in enumerate(self.current_candidates):
            self.current_A[:, count] = self.observe(loc)
        self.current_A[:, len(self.current_candidates)] = 1

    def set_check_mark(self):
        self.hist_candidates.append(self.current_candidates)
        self.hist_candidates_intervals.append(self.current_candidates_intervals)
        self.hist_PRIS_ItN.append(self.current_PRIS_ItN)
        self.hist_check_mark_id.append(self.current_check_mark_id)
        self.current_check_mark_id += 1

    def show_attributes(self):
        '''
        This is for a convenient check of the attributes.
        :return:  display a 000000list of attributes with the corresponding values, except for the long ones.
        '''
        for key, value in self.__dict__.items():
            if key is "current_candidates":
                print(key + ":  [hidden];")
            elif key is "hist_candidates":
                print(key + ":  [hidden]")
            elif key is "current_A":
                print(key + ":  [hidden];")
            else:
                print(key + ":  " + str(value))


class LinBreg:
    import time
    def __init__(self):
        self.id = []  # ID unique to PRIS object
        # solve for x from Ax = b.
        self.A = 0  # sensing matrix.
        self.x = 0
        self.b = 0  # observation vector.
        self.flag_stop = False  # flag to stop optimization iteration.
        self.maxit = 2000  # maximum iteration steps.
        self.debug_it_int = 1
        self.flag_positivity = True
        self.it_check_rem = 1
        self.iterations = list()
        self.hist_res = list()
        self.hist_resDrop = list()
        self.save_obj_int = 100
        self.bg = list()
        self.alpha = 1
        self.debug = False
        self.deep_debug = False
        self.save = True
        self.obs_dim0 = 0
        self.obs_dim1 = 0

        self.kick = self.Kick(self)
        
        self.A_dir = ''  # directory to store sensing matrix when saving
        try:
            with open("../../PyPRIS_Scratch/saved_objects/PyPRIS_{}_SensingMx.file".format(self.id), "wb") as f:
                pickle.dump(self.A, f, pickle.HIGHEST_PROTOCOL)
            except OSError:
                print ("Failed to write sensing matrix to directory %s " % path_s)
            else:
                print ("Successfully wrote sensing matrix to directory %s " % path_s)

    class Kick:
        def __init__(self, LinBreg):
            self.parent = LinBreg
            self.flag = False
            self.ints = 10  # number of iterations between kicking evaluation.
            self.reference = copy.deepcopy(self.parent.x)
            self.option = False
            self.thres = 1e-10
            self.refnorm = np.max(self.parent.x)
            self.hist_refnorm = list()
            self.hist_eval_counts = list()

        def evaluation(self, it_count):
            self.refnorm = np.linalg.norm(self.parent.x - self.reference)
            if self.refnorm < self.thres:
                # flip the kicking flag to "True" with positive evaluation
                self.flag = True
            else:
                # flip the kicking flag to "False" with negative evaluation.
                self.flag = False

            self.hist_refnorm.append(self.refnorm)
            self.hist_eval_counts.append([it_count, self.flag])

        def go(self):
            # execute kicking
            # kicking execution only modifies the domains of the step size 
            # and split domains in a binary way.
            i0 = np.where(self.parent.x == 0)  # zero entries on x [x entries where there kicking is in demand]
            i1 = np.where(
                self.parent.x != 0)  # none zero entries in x. [x entries where there is a value, no need for kicking]
            si = (self.parent.mu * np.sign(self.parent.erpj[i0]) - self.parent.cumerr[i0]) / self.parent.erpj[
                i0]  # stepsie for entries that needs kicking
            self.parent.stepsize[i0] = np.min(si)
            self.parent.stepsize[i1] = 1
            # reset kick.flag to False and wait for the flip 
            # from the next positive kicking evaluation
            self.flag = False
            # update kick.reference for follow-up kicking evaluation
            self.reference = copy.deepcopy(self.parent.x)

    def get_ready(self):
        self.x = np.zeros(self.A.shape[1])
        self.stepsize = np.ones(self.x.shape)  # stepsize.
        self.er = np.zeros(self.b.shape)
        self.erpj = np.zeros(self.x.shape)
        self.cumerr = np.zeros(self.x.shape)
        self.recb = np.zeros(self.b.shape)
        
        import os
        # define the name of the directory to be created.
        path_0 = "../../PyPRIS_Scratch/"
        path_d = "../../PyPRIS_Scratch/debug_output"
        path_s = "../../PyPRIS_Scratch/saved_objects"
        try:
            if not os.path.exists(path_0):
                os.mkdir(path_0)
        except OSError:
            print ("Creation of the directory %s failed" % path_0)
        else:
            print ("Successfully created Scratch directory %s " % path_0)

        if self.save is True:
            try:
                if not os.path.exists(path_s):
                    os.mkdir(path_s)

            except OSError:
                print ("Creation of the directory %s failed" % path_s)
            else:
                print ("Successfully created Object-saving directory %s " % path_s)

        if self.debug is True:
            try:
                if not os.path.exists(path_d):
                    os.mkdir(path_d)
            except OSError:
                print ("Creation of the directory %s failed" % path_d)
            else:
                print ("Successfully created Debug directory %s " % path_d)

    def shrink(self, sk):
        sk[np.where((sk >= -self.mu) * (sk <= self.mu))] = 0
        sk[np.where(sk > self.mu)] -= self.mu
        sk[np.where(sk < -self.mu)] += self.mu
        return sk

    def go(self):
        t1 = time.time()
        it_count = 0
        self.hist_res.append(0)
        self.hist_resDrop.append(0)
        self.iterations.append(it_count)
        self.bg.append(self.x[self.x.size - 1])
        # main linearized bregman iteration with kicking option.
        while self.flag_stop is False:
            # incrementation of the iteration number.
            it_count += 1

            # calculate distance (error)
            self.recb = np.dot(self.A, self.x)
            self.er = self.b - self.recb
            if self.deep_debug is True: self.debug_output(it_count, appstr='_a_er_updated')

            # perform back projection of the error ('adding the errors back').
            self.erpj = np.dot(self.er, self.A)
            if self.deep_debug is True: self.debug_output(it_count, appstr='_b_erpj_updated')

            # check if kicking is needed 
            #
            # "Kicking" rescales the backprojected error (self.erpj) with two different stepsizes
            # we'll have stepsize > 1 for kicking area, and stepsize = 1 for non-kicking area. 
            # kicking boosts the tip of the cumulated backprojected errors towards the shrinkage 
            # threshold.
            # In this implementation, the effect of kicking is realized throught a modified. 
            # distribution of stepsizes (self.stepsize). 
            self.stepsize = np.ones(self.x.shape)  # [Note: this step involves some redundancy]
            if np.remainder(it_count, self.kick.ints) == 0:
                self.kick.evaluation(it_count)
                if self.deep_debug is True: self.debug_output(it_count, appstr='_c1_kicking_evaluated')
                # kick if we get a positive kicking ealuation.
                if self.kick.flag is True:
                    self.kick.go()
                if self.deep_debug is True: self.debug_output(it_count, appstr='_c2_kicking_updated')

            # get the acumulation of the back projected error.
            self.cumerr += self.erpj * self.stepsize
            if self.deep_debug is True: self.debug_output(it_count, appstr='_d_cumerr_updated')

            # perform positivity constraint:
            if self.flag_positivity is True: self.x[np.where(self.x < 0)] = 0
            if self.deep_debug is True: self.debug_output(it_count, appstr='_e_positivity_updated')

            # shrinkage to update the candidate coefficients.
            self.x = copy.deepcopy(self.cumerr)
            if self.deep_debug is True: self.debug_output(it_count, appstr='_f_x_copied')
            self.x = self.alpha * self.shrink(self.x)
            if self.deep_debug is True: self.debug_output(it_count, appstr='_g_x_updated')

            # decide on the termination of iterations.
            if it_count > self.maxit: self.flag_stop = True

            # update the quantities for status tracking purposes.
            self.track_status(it_count, self.er)

            # check intermediate outputs. (Valid under debug mode).
            self.debug_output(it_count, appstr='_h_track_status_updated')

            # save object into separate file every assigned step
            self.save_obj(it_count, self.save_obj_int)

    def save_obj(self, currit, step):
        if self.save is True:
            if currit % step == 1:
                self.A = 0
                with open("../../PyPRIS_Scratch/saved_objects/PyPRIS_{}_{}.file".format(self.id, currit), "wb") as f:
                    pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
                    print ("Successfully saved Linbreg ID {} at iteration {} to directory.".format(self.id, currit))
                with open('../../PyPRIS_Scratch/saved_objects/PyPRIS_{}_SensingMx.file'.format(self.id), "rb") as s:
                    self.A = pickle.load(s)
                    
    def candidate_vis(self):
        intervals = self.candidate_intervals
        locs = list(zip(*self.candidate_coords))
        dims = list()
        minimals = list()
        maximums = list()

        for inds, interval in zip(locs, intervals):
            maximum = np.max(inds)
            minimal = np.min(inds)
            dims.append(1 + np.int((maximum - minimal) // interval))
            minimals.append(minimal)
            maximums.append(maximum)

        vis = np.zeros(dims)
        for coords, intensity in zip(self.candidate_coords, self.x[0:len(self.x) - 1]):
            vis[coords[0] - minimals[0] - 1, int((coords[1] - minimals[1]) // intervals[1]), int(
                (coords[2] - minimals[2]) // intervals[2])] = intensity

        return vis

    def debug_output(self, it_count, appstr):
        # Generate intermediate output under debug mode.
        if self.debug is True:
            if np.remainder(it_count, self.debug_it_int) == self.it_check_rem:
                print('intermediate output it#' + str(it_count))
                vis = self.candidate_vis()
                nrow = 4
                ncol = 5
                plt.figure(figsize=(13, 8))
                plt.subplot(nrow, ncol, 1)
                if len(vis) is 0:
                    t = plt.title("No signal yet.")
                else:
                    plt.imshow(np.mean(vis, axis=0))
                    t = plt.title('XY-plane projection')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 2)
                plt.imshow(self.b.reshape(self.obs_dim0, self.obs_dim1))
                t = plt.title('input blur')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 3)
                plt.imshow(self.recb.reshape(self.obs_dim0, self.obs_dim1))
                t = plt.title('recovered blur')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 4)
                plt.plot(self.recb.ravel(), '.')
                t = plt.title('recovered obsrvation plot')
                t.set_position([.5, 1.15])

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
                plt.plot([0, len(self.cumerr.ravel())], [self.mu, self.mu], 'r')
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
                plt.plot(self.stepsize.ravel(), '.')
                t = plt.title('stepsize distribution')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 12)
                plt.plot(self.iterations, np.log(self.hist_res), '.')
                t = plt.title('Log(L2(res))')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 13)
                plt.plot(self.er.ravel(), '.')
                t = plt.title('err')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 14)
                if len(self.kick.hist_eval_counts) > 2:
                    t = list(zip(*self.kick.hist_eval_counts))
                    plt.scatter(t[0], t[1], c=t[1])
                    t = plt.title('kicking history')
                    t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 16)
                plt.text(0, 0.9, 'mu: ' + str(np.floor(self.mu)), fontsize=16)
                plt.text(0, 0.6, 'stepsize: ' + str(np.floor(self.stepsize)), fontsize=16)
                plt.text(0, 0.3, 'current kicking flag is: ' + str(self.kick.flag), fontsize=16)
                plt.text(0, 0, 'current figure: plots_it' + str(it_count) + appstr, fontsize=16)
                plt.axis('off')

                plt.subplots_adjust(top=0.95, left=0, right=1, bottom=0, wspace=0.5, hspace=1)
                plt.savefig(
                    '../../PyPRIS_Scratch/debug_output/PyPRIS_{}_plots_it{}{}.png'.format(self.id, it_count, appstr),
                    dpi=300, figsize=(100, 80))
                plt.close()

    def track_status(self, it_count, er):
        self.hist_res.append(np.linalg.norm(er))
        self.hist_resDrop.append((self.hist_res[it_count] - self.hist_res[it_count - 1]) / self.hist_res[it_count - 1])
        self.iterations.append(it_count)
        self.bg.append(self.x[self.x.size - 1])


def loadPyPRIS(PRIS_iter, step):
    with open('../../PyPRIS_Scratch/saved_objects/PyPRIS_{}_{}.file'.format(PRIS_iter, step), "rb") as f:
        PyPRIS = pickle.load(f)
    with open('../../PyPRIS_Scratch/saved_objects/PyPRIS_{}_SensingMx.file'.format(PRIS_iter), "rb") as s:
        PyPRIS.A = pickle.load(s)
    return PyPRIS
