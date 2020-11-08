import numpy as np
import time
import copy
import pickle
import joblib
import matplotlib

from PyPRIS.utils import *
from PyPRIS.candidate_screening import *


try:
    from matplotlib import pyplot as plt
    plt.switch_backend('agg')
except RuntimeError:
    pass

#  Authors: Xiyu Yi, Xingjia Wang @ UCLA, 2019.
#  PI: Shimon Weiss, Department of Chemistry and Biochemistry, UCLA.
#  Xiyu continued to work on the codes after July 8th of 2020,
#  when she started to become an employee at Lawrence Livermore National Laboratory.


class PyPRIS:
    
    def __init__(self):
        self.name = 'PyPRIS object'
        self.positivity = True  # positivity constraint.
        self.observation = np.ndarray(0)  # this should be the observation vector
        self.current_relReF = list([1, 2, 2])    # current relative refinement factor. This is the relative refinement
        #                               to be or have been performed
        #                               for this round of pris as compared to the last round of pris.
        self.current_PRIS_ItN = 0  # current PRIS iteration count, starting from 0.
        self.current_A = np.ndarray(0)  # current sensing matrix.
        self.current_candidates = list()  # current pool of candidates.
        self.current_candidates_intervals = list()  # current intervals of neighboring candidate voxels.
        self.current_check_mark_id = 0

        self.hist_candidates = list()  # keep a record of the full history.
        self.hist_candidates_intervals = list()  # keep a record of the full history.
        self.hist_PRIS_ItN = list()
        self.hist_check_mark_id = list()  # checkmark ID. ascending each time after you make a check_mark
        
        self.ifsave = True
        self.path_s = "./saved_objects"
        self.expansion = False
        self.species_n = 1
        self.top_candidates = False
        self.top_candidates_N = 500

    def save(self):
        import os
        try:
            if not os.path.exists(self.path_s):
                os.mkdir(self.path_s)
        except OSError:
            print ("Creation of the directory %s failed" % self.path_s)
        else:
            print ("Successfully created Scratch directory %s " % self.path_s)
        if self.ifsave is True:
            self.current_A = np.ndarray(0)
            with open("{}/PyPRIS_pris{}.file".format(self.path_s, self.current_PRIS_ItN), "wb") as f:
                pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
                print ("Successfully saved PyPRIS ID {} to directory.".format(self.current_PRIS_ItN))
                    
    def prep_for_new_refinement(self):
        self.hist_candidates.append(copy.deepcopy(self.current_candidates))
        self.hist_candidates_intervals.append(copy.deepcopy(self.current_candidates_intervals))
        self.hist_PRIS_ItN.append(copy.deepcopy(self.current_PRIS_ItN))
        self.set_check_mark()

    def refine_candidates(self, cs_solver):
    # this will take the current candidates and the result in the compressive sensing solver object,
    # and generate refined pool of candidates.
        # create a check mark.
        self.set_check_mark()
        # Get the non_zero_coordinates from the existing cs_solver results.
        non_zero_inds = np.argwhere(cs_solver.x[0:-1] > 0)
        if self.top_candidates is True:
            print("refine candidates with only top candidates")
            maxN = np.min([self.top_candidates_N, len(non_zero_inds)])
            non_zero_inds = cs_solver.x[0:-1].argsort()[-maxN:][::-1]

        non_zero_coordinates = [self.current_candidates[i] for i in list(non_zero_inds.ravel())]
        self.current_candidates_intervals = copy.deepcopy(
            [pre / ref for pre, ref in zip(self.current_candidates_intervals, self.current_relReF)]
        )
        current_interval = self.current_candidates_intervals
    # get new coordinates with 2-fold refinement.
        new_coords = list()
        print('expansion is ' + str(self.expansion))
        for i in non_zero_coordinates:
            if self.expansion is False:
                if self.species_n == 2:
                    extra_coords = [[i[0], i[1] - current_interval[1] / 2, i[2] - current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] - current_interval[1] / 2, i[2] + current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] + current_interval[1] / 2, i[2] - current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] + current_interval[1] / 2, i[2] + current_interval[2] / 2, i[3]]]
                else:
                    extra_coords = [[i[0], i[1] - current_interval[1] / 2, i[2] - current_interval[2] / 2], \
                                    [i[0], i[1] - current_interval[1] / 2, i[2] + current_interval[2] / 2], \
                                    [i[0], i[1] + current_interval[1] / 2, i[2] - current_interval[2] / 2], \
                                    [i[0], i[1] + current_interval[1] / 2, i[2] + current_interval[2] / 2]]
            else:
                if self.species_n == 2:
                    extra_coords = [[i[0], i[1] -     current_interval[1] / 2, i[2] -     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] -     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] -     current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2, i[3]], \
                                    \
                                    [i[0], i[1] -     current_interval[1] / 2, i[2] +     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] +     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] -     current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2, i[3]], \
                                    \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] -     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] -     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2, i[3]], \
                                    \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] +     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] +     current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2, i[3]], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2, i[3]]]
                else:
                    extra_coords = [[i[0], i[1] -     current_interval[1] / 2, i[2] -     current_interval[2] / 2], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] -     current_interval[2] / 2], \
                                    [i[0], i[1] -     current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2], \
                                    \
                                    [i[0], i[1] -     current_interval[1] / 2, i[2] +     current_interval[2] / 2], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] +     current_interval[2] / 2], \
                                    [i[0], i[1] -     current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2], \
                                    [i[0], i[1] - 3 * current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2], \
                                    \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] -     current_interval[2] / 2], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] -     current_interval[2] / 2], \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] - 3 * current_interval[2] / 2], \
                                    \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] +     current_interval[2] / 2], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] +     current_interval[2] / 2], \
                                    [i[0], i[1] +     current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2], \
                                    [i[0], i[1] + 3 * current_interval[1] / 2, i[2] + 3 * current_interval[2] / 2]]

            for i1 in extra_coords:
                if i1 not in new_coords:
                    new_coords.append(i1)

        self.current_candidates = copy.deepcopy(new_coords)

        # set a check mark for tracking purposes.
        self.set_check_mark()
        
    def candidate_pop(self, cs_solver, thres=0):
        """
        pop out candidates with amplitudes below or equal to the threshold.
        """
        # Note: 
        # maitain the separation of pypris from the cs_solver type. 
        # don't modify cs_solver in pypris methods. 
        survival_inds = np.argwhere(cs_solver.x[0:len(cs_solver.x) - 1] > thres)
        if len(survival_inds)>0:
            survival_coordinates = [self.current_candidates[i] for i in list(survival_inds.ravel())]
            self.current_candidates = copy.deepcopy(survival_coordinates)
            self.survival_inds = survival_inds

    def generate_sensing_mx(self):
        print("----------- Generate sensing matrix:")
        print("            Matrix size:",str(len(self.observation)),' observation pixels ')
        print("                        ",str(len(self.current_candidates)),' candidates ')
        self.current_A = np.ndarray([len(self.observation), len(self.current_candidates) + 1])
        for count, loc in enumerate(self.current_candidates):
            if self.species_n == 2:
                self.current_A[:, count] = self.observe(loc[0:3],loc[3])
            else:
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

    def prep_candit_rolling_pop(self):
        """
        this function takes a pypris object in, and returns a modified version.
        :return:
        """
        self = utils.candidates_init_rolling_pop(self)
        return


class LinBreg:
    def __init__(self, PyPRIS_n):
        import time
        self.auto_mu = False
        self.auto_bg = False
        self.auto_mu_fold = 1
        self.auto_bg_fold = 2
        self.PyPRIS_iter = []  # Associated PyPRIS iter number
        self.PyPRIS_name = PyPRIS_n # Associated PyPRIS name
        self.path_0 = "../PyPRIS_Scratch"
        # solve for x from Ax = b.
        self.A = 0  # sensing matrix.
        self.x = np.zeros(0)
        self.b = 0  # observation vector.
        self.flag_stop = False  # flag to stop optimization iteration.
        self.maxit = 2000  # maximum iteration steps.
        self.debug_it_int = 1
        self.flag_positivity = True
        self.flag_bg_allow_negative = True
        self.it_check_rem = 1
        self.iterations = list()
        self.hist_res = list()
        self.hist_delta_res = list()
        self.hist_percent_delta_res = list()
        self.stopping_moni_start = 1000
        self.it_count = 0
        self.candidate_coords = list()
        # now initialize a threshold value for abs(log(abs(x)))*sign(x) with x=hist_delta_res; it suppose to be a negative value.
        self.stopping_loghistpercdelres_thres = -11; # iteration will stop with the value is below this value.
        self.stopping_loghistpercdelres = np.Inf;
        self.save_obj_int = 100
        self.bg = list()
        self.alpha = 1
        self.save = True
        self.obs_dim0 = 0
        self.obs_dim1 = 0
        self.kick = self.Kick(self)
        self.A_dir = ''  # directory to store sensing matrix when saving
        self.mu_reference = 0


    class Kick:
        def __init__(self, LinBreg):
            self.parent = LinBreg
            self.flag = False
            self.ints = 10  # number of iterations between kicking evaluation.
            self.option = False
            self.thres = 1e-10
            tm = list(self.parent.x.ravel())
            tm.append(-np.Inf)
            self.refnorm = np.max(tm)
            self.hist_refnorm = list()
            self.hist_eval_counts = list()

        def set_reference(self):
            self.reference = copy.deepcopy(self.parent.x)

        def evaluation(self, it_count):
            self.refnorm = np.linalg.norm(self.parent.x - self.reference)
            if self.refnorm < self.thres:
                # flip the kicking flag to "True" with positive evaluation
                print("----- kick evaluation:")
                print("----- distance of x change from "+str(self.ints)+" steps before is:"+str(self.refnorm))
                print("----- Threshold is:"+str(self.thres))
                print("------------------------------- Kick!")
                self.flag = True
            else:
                # flip the kicking flag to "False" with negative evaluation.
                print("----- kick evaluation:")
                print("----- distance of x change from "+str(self.ints)+" steps before is:"+str(self.refnorm))
                print("----- Threshold is:"+str(self.thres))
                print("---------------- No kick.")
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
            si = (self.parent.mu * np.sign(self.parent.respj[i0]) - self.parent.cumres[i0]) / self.parent.respj[
                i0]  # stepsie for entries that needs kicking
            if len(si) > 0:
                self.parent.stepsize[i0] = np.min(si)
            self.parent.stepsize[i1] = 1
            # reset kick.flag to False and wait for the flip 
            # from the next positive kicking evaluation
            self.flag = False
            # update kick.reference for follow-up kicking evaluation
            self.reference = copy.deepcopy(self.parent.x)

    def apply_mask(self, mask):
    # this will take the current candidates and the result in the compressive sensing solver object,
    # and generate  pool of candidates.        
        vis1 = self.candidate_vis()

        # when the mask value is 1, don't change the value. when the mask value is 0, set it to 0.
        vis1 = mask * vis1
        self.x = copy.deepcopy(self.candidate_vis_inv(vis1))

    def match_popped_candidates(self, survival_inds):
        xtemp = np.ndarray(len(survival_inds) + 1)
        for i, iv in enumerate(survival_inds):
            xtemp[i] = self.x[iv]
        xtemp[-1] = self.x[-1]
        self.x = xtemp

    def candidate_pool_thinning(self, opts):
        if isinstance(opts, ProjFilter2D_Opts):
            v = self.candidate_vis()
            prj = copy.deepcopy(np.mean(v, axis=0))
            ub = np.max(prj.ravel())
            lb = np.min(prj.ravel())
            thres = lb + opts.relative_lower_bound_2Dproj*(ub - lb)
            prj[np.where(prj < thres)] = 0
            prj[np.where(prj >= thres)] = 1
            for i in np.arange(0,v.shape[0]):
                v[i,:,:] = copy.deepcopy(prj)
            mask = v
            self.mask_2D = prj

        if isinstance(opts, RelativeValueFilter_Opts):
            v = self.candidate_vis()
            ub = np.max(prj.ravel())
            lb = np.min(prj.ravel())
            thres = lb + opts.relative_lower_bound_FullPool * (ub - lb)
            v[np.where(v < thres)] = 0
            v[np.where(v >= thres)] = 1
            mask = v

        return mask

    def get_ready(self):
        import os
        print("linbreg getready:")
        self.it_count = -1
        self.x = np.zeros(self.A.shape[1])
        self.stepsize = np.ones(self.x.shape)  # stepsize.
        self.res = np.zeros(self.b.shape)
        self.respj = np.zeros(self.x.shape)
        self.cumres = np.zeros(self.x.shape)
        self.recb = np.zeros(self.b.shape)
        self.path_s = self.path_0 + "/saved_objects"
        self.path_d = self.path_0 + "/debug_output"
        self.kick.set_reference() # set a kick reference
        print('stopping threshold: '+str(self.stopping_loghistpercdelres_thres))
        print('alpha: '+str(self.alpha))
        # calculate the initial back projected errors to determine mu
        # define the name of the directory to be created.
        import os
        try:
            if not os.path.exists(self.path_0):
                os.mkdir(self.path_0)
        except OSError:
            print ("Creation of the directory %s failed" % self.path_0)
        else:
            print ("Successfully created Scratch directory %s " % self.path_0)

        if self.save is True:   
            try:
                if not os.path.exists(self.path_s):
                    os.mkdir(self.path_s)

                try:
                    with open("{}/PyPRIS_{}_{}_SensingMx.file".format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter), "wb") as f:
                        joblib.dump(self.A, f, pickle.HIGHEST_PROTOCOL)
                except OSError:
                    print ("Failed to write sensing matrix to directory %s " % self.path_s)
                else: 
                    print ("Successfully wrote sensing matrix to directory %s " % self.path_s)

                self.A = 0 # remove sensing matrix because it is too big for pickle
                with open("{}/PyPRIS_{}_{}_{}.file".format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter, 0), "wb") as f:
                    pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
                    print ("Successfully saved Linbreg ID {} at iteration {} to directory.".format(self.PyPRIS_iter, 0))
                with open('{}/PyPRIS_{}_{}_SensingMx.file'.format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter), "rb") as s:
                    self.A = joblib.load(s)  # load sensing matrix back

            except OSError:
                print ("Creation of the directory %s failed" % self.path_s)
            else:
                print ("Successfully created Object-saving directory %s " % self.path_s)

        if self.debug is True:
            try:
                if not os.path.exists(self.path_d):
                    os.mkdir(self.path_d)
            except OSError:
                print ("Creation of the directory %s failed" % self.path_d)
            else:
                print ("Successfully created Debug directory %s " % self.path_d)

        if self.auto_bg is True:
            print('setting auto background component transfer function')
            pp = np.max(np.dot(self.b, self.A[:, :-1]))
            bb = np.max(np.dot(self.b, self.A[:, -1]))
            r = pp/bb
            self.A[:, -1] = self.A[:, -1] * r * self.auto_bg_fold

            # unfinished...

        if self.auto_mu is True:
            print('Setting AutoMu')
            pp = np.dot(self.b, self.A)
            self.mu = np.max(pp.ravel())*self.auto_mu_fold

    def shrink(self, sk):
        sk[np.where((sk >= -self.mu) * (sk <= self.mu))] = 0
        sk[np.where(sk > self.mu)] -= self.mu
        sk[np.where(sk < -self.mu)] += self.mu
        return sk

    def go(self):
        t1 = time.time()
        self.hist_res.append(0)
        self.hist_delta_res.append(0)
        self.hist_percent_delta_res.append(0)
        self.iterations.append(self.it_count)
        self.bg.append(self.x[self.x.size - 1])
        # main linearized bregman iteration with kicking option.
        while self.flag_stop is False:
            # incrementation of the iteration number.
            self.it_count += 1
            it_count = self.it_count

            # calculate distance (residuals)
            self.recb = np.dot(self.A, self.x)
            self.res = self.b - self.recb
            if self.deep_debug is True: self.debug_output(it_count, appstr='_a_res_updated')

            # perform back projection of the residuals ('adding the residuals back').
            self.respj = np.dot(self.res, self.A)
            if self.deep_debug is True: self.debug_output(it_count, appstr='_b_respj_updated')

            # check if kicking is needed
            #
            # "Kicking" rescales the back projected residual (self.respj) with two different stepsizes
            # we'll have stepsize > 1 for kicking area, and stepsize = 1 for non-kicking area.
            # kicking boosts the tip of the cumulated backprojected residuals towards the shrinkage
            # threshold.
            # In this implementation, the effect of kicking is realized throught a modified.
            # distribution of stepsizes (self.stepsize).
            self.stepsize = np.ones(self.x.shape)  # [Note: this step involves some redundancy]
            if np.remainder(it_count, self.kick.ints) == 1:
                self.kick.evaluation(it_count)
                self.kick.set_reference()
                if self.deep_debug is True: self.debug_output(it_count, appstr='_c1_kicking_evaluated')
                # kick if we get a positive kicking ealuation.
                if self.kick.flag is True:
                    print("Now kick.")
                    self.kick.go()
                    self.kick.set_reference() # update the kick reference
                if self.deep_debug is True: self.debug_output(it_count, appstr='_c2_kicking_updated')

            # get the acumulation of the back projected residuals.
            self.cumres += self.respj * self.stepsize
            if self.deep_debug is True: self.debug_output(it_count, appstr='_d_cumres_updated')

            # shrinkage to update the candidate coefficients.
            self.x = copy.deepcopy(self.cumres)
            if self.deep_debug is True: self.debug_output(it_count, appstr='_f_x_copied')
            self.x = self.alpha * self.shrink(self.x)
            if self.deep_debug is True: self.debug_output(it_count, appstr='_g_x_updated')

            # perform positivity constraint:
            bg = self.x[-1]
            if self.flag_positivity is True:
                self.x[np.where(self.x < 0)] = 0
                if self.flag_bg_allow_negative is True:
                    self.x[-1]=bg # release the background component from the positivity constraint

            if self.deep_debug is True: self.debug_output(it_count, appstr='_e_positivity_updated')

            # update the quantities for status tracking purposes.
            self.track_status(it_count, self.res)

            # decide on the termination of iterations.
            # set termination signal if maximum iteration is reached:
            if it_count >= self.maxit: self.flag_stop = True
            # set termination signal if stopping criteria is met:
            if self.stopping_loghistpercdelres < self.stopping_loghistpercdelres_thres and it_count > self.stopping_moni_start: 
                self.flag_stop = True
                print('stopping criteria fulfilled')
            # check intermediate outputs. (Valid under debug mode).
            self.debug_output(it_count, appstr='_h_track_status_updated')
            # save object into separate file every assigned step
            self.save_obj(it_count, self.save_obj_int)
                
    def save_obj(self, currit, step):
        import sys
        if self.save is True:
            if currit % step == 1:
                print("current iteration " + str(self.it_count))
                print("now start saving objs")
                self.A = 0 # remove sensing matrix because it is too big.
                setattr(sys.modules[__name__], 'Kick', self.Kick)
                with open("{}/PyPRIS_{}_{}_{}.file".format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter, currit), "wb") as f:
                    pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
                    print ("Successfully saved Linbreg ID {} at iteration {} to directory.".format(self.PyPRIS_iter, currit))
                with open('{}/PyPRIS_{}_{}_SensingMx.file'.format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter), "rb") as s:
                    self.A = joblib.load(s) # load sensing matrix back
            elif self.flag_stop is True:
                self.A = 0 # remove sensing matrix because it is too big
                with open("{}/PyPRIS_{}_{}_{}.file".format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter, currit), "wb") as f:
                    pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
                    print ("Successfully saved Linbreg ID {} at iteration {} to directory.".format(self.PyPRIS_iter, currit))
                with open('{}/PyPRIS_{}_{}_SensingMx.file'.format(self.path_s, self.PyPRIS_name, self.PyPRIS_iter), "rb") as s:
                    self.A = joblib.load(s) # load sensing matrix back

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
            tp = vis[int((coords[0] - minimals[0] - 1) // intervals[0]), int((coords[1] - minimals[1]) // intervals[1]), int(
                (coords[2] - minimals[2]) // intervals[2])]
            vis[int((coords[0] - minimals[0] - 1) // intervals[0]), int((coords[1] - minimals[1]) // intervals[1]), int(
                (coords[2] - minimals[2]) // intervals[2])] = np.max([intensity, tp])

        return vis
    
    def candidate_vis_inv(self, vis):
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
        
        new_x = []
        for coords, intensity in zip(self.candidate_coords, self.x[0:len(self.x) - 1]):
            new_x.append(vis[coords[0] - minimals[0] - 1, int((coords[1] - minimals[1]) // intervals[1]), int(
                (coords[2] - minimals[2]) // intervals[2])] )
                    
        new_x.append(self.x[-1])
        
        return np.array(new_x)

    def debug_output(self, it_count, appstr):
        # Generate intermediate output under debug mode.
        if self.debug is True:
            if np.remainder(it_count, self.debug_it_int) == self.it_check_rem or self.flag_stop is True:
                print('intermediate output it#' + str(it_count))
                vis = self.candidate_vis()
                nrow = 5
                ncol = 5
                plt.figure(figsize=(14, 9))
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
                t = plt.title('recovered observation plot')
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
                plt.plot(self.cumres.ravel(), '.')
                t = plt.title('cum res')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 8)
                plt.plot(self.cumres.ravel(), '.')
                plt.plot([0, len(self.cumres.ravel())], [self.mu, self.mu], 'r')
                t = plt.title('cum res and mu')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 9)
                plt.plot(self.respj.ravel(), '.')
                t = plt.title('residual back projection (respj)')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 10)
                plt.plot(self.b.ravel(), '.')
                t = plt.title('input observation plot')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 11)
                plt.plot(self.iterations, (self.hist_res), '.')
                t = plt.title('L2(res)')
                t.set_position([.5, 1.15])

                plt.subplot(nrow, ncol, 12)
                plt.plot(self.iterations, self.hist_delta_res, '.')
                t = plt.title('delta L2(res)')
                t.set_position([.5, 1])
                
                plt.subplot(nrow, ncol, 13)
                plt.plot(self.iterations, self.hist_percent_delta_res, '.')
                t = plt.title('percent delta L2(res)')
                t.set_position([.5, 1])
                
                plt.subplot(nrow, ncol, 14)
                plt.plot(self.res.ravel(), '.')
                t = plt.title('res')
                t.set_position([.5, 1])

                plt.subplot(nrow, ncol, 15)
                plt.plot(self.stepsize.ravel(), '.')
                t = plt.title('stepsize distribution')
                t.set_position([.5, 1])
                

                plt.subplot(nrow, ncol, 16)
                plt.plot(self.iterations, abs(np.log(abs(np.asarray(self.hist_res)))) \
                                          *np.sign(self.hist_res), \
                         '.')
                t = plt.title('abs(log(abs(x)))*sign(x) \n x = Log(L2(res))')
                t.set_position([.5, 1])

                plt.subplot(nrow, ncol, 17)
                plt.plot(self.iterations, abs(np.log(abs(np.asarray(self.hist_delta_res)))) \
                                          *np.sign(self.hist_delta_res), \
                         '.')
                t = plt.title('abs(log(abs(x)))*sign(x) \n x = delta_L2(res)')
                t.set_position([.5, 1])

                plt.subplot(nrow, ncol, 18)
                plt.plot(self.iterations, abs(np.log(abs(np.asarray(self.hist_percent_delta_res)))) \
                                          *np.sign(self.hist_percent_delta_res), \
                         '.')
                plt.plot([self.iterations[0],self.iterations[-1]], \
                         [self.stopping_loghistpercdelres_thres,self.stopping_loghistpercdelres_thres],'r')
                t = plt.title('abs(log(abs(x)))*sign(x) \n x = percent_delta_L2(res)')
                t.set_position([.5, 1])

                plt.subplot(nrow, ncol, 19)
                try:
                    plt.hist(self.res.ravel(), 100)
                    t = plt.title('histogram of residual')
                except:
                    t = plt.title('histogram of residual not available')
                    pass
                t.set_position([.5, 0.5])


                plt.subplot(nrow, ncol, 20)
                if len(self.kick.hist_eval_counts) > 2:
                    t = list(zip(*self.kick.hist_eval_counts))
                    plt.scatter(t[0], t[1], c=t[1])
                    t = plt.title('kicking history')
                    t.set_position([.5, 1])      
                
                plt.subplot(nrow, ncol, 24)
                plt.text(0, 1, 'mu: ' + str(np.floor(self.mu)), fontsize=14)
                plt.text(0, 0.8, 'stepsize: ' + str(np.floor(self.stepsize)), fontsize=14)
                plt.text(0, 0.6, 'current kicking flag is: ' + str(self.kick.flag), fontsize=14)
                plt.text(0, 0.4, 'current figure:', fontsize=14)
                plt.text(0, 0.2, 'plots_it' + str(it_count) + appstr, fontsize=14)
                plt.axis('off')
                plt.subplots_adjust(top=0.95, left=0.1, right=0.9, bottom=0.1, wspace=0.5, hspace=1)
                plt.savefig(
                    '{}/PyPRIS_{}_{}_plots_it{}{}.png'.format(self.path_d, self.PyPRIS_name, self.PyPRIS_iter, it_count, appstr),
                    dpi=300, figsize=(100, 80))
                plt.close()

    def track_status(self, it_count, res):
        self.hist_res.append(np.linalg.norm(res))
        self.hist_delta_res.append((self.hist_res[it_count] - self.hist_res[it_count - 1]))
        self.hist_percent_delta_res.append( \
                                           (self.hist_res[it_count] - self.hist_res[it_count - 1]) \
                                           /self.hist_res[it_count - 1] \
                                          )
        self.iterations.append(it_count)
        self.bg.append(self.x[self.x.size - 1])
        stopping_tag = copy.deepcopy(self.hist_percent_delta_res[-1])
        self.stopping_loghistpercdelres = abs(np.log(abs(stopping_tag)))*np.sign(stopping_tag)


def loadCSSolver(path, PyPRIS_name, PyPRIS_SensMx_name):
    with open('{}/{}.file'.format(path, PyPRIS_name), "rb") as f:
        linbreg = pickle.load(f) #the loaded object is a LinBreg object

    try:
        with open('{}/{}.file'.format(path, PyPRIS_SensMx_name), "rb") as s:
            linbreg.A = joblib.load(s)
    except:
        print("the following file not found:")
        print('{}/{}.file'.format(path, PyPRIS_SensMx_name))
        print("No sensing matrix available, the relevant visualization will be omitted.")
        pass
    return linbreg


def loadPyPRIS(path, PyPRIS_name):
    with open('{}/{}.file'.format(path, PyPRIS_name), "rb") as f:
        pris = pickle.load(f) #the loaded object is a PyPRIS object
#     with open('{}/{}.file'.format(path, PyPRIS_SensMx_name), "rb") as s:
#         pris.current_A = joblib.load(s)
    return pris

def get_ticket(ticket_path):
    with open(ticket_path, "rb") as f:
        ticket = pickle.load(f)
    return ticket
