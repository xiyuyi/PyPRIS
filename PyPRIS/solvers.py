import numpy as np
import torch


class PRIS:
    """

    """
    def __init__(self):
        self.type = 'PRIS solver'
        self.refinement_number = 5  # default total number of refinement
        self.refinement_factor_list = [1, 2, 2, 2, 2]  # default list of refinement factors
        self.lasso_solver = 'lbreg'  # default choice of lasso solver. Only available one currently.
        self.dual_channel = 'Yes'  # default setting of the number of observation channels.
        self.A = []
        self.b = []
        self.x = torch.empty(1, 1, dtype=torch.float)
        self.lbreg_opts = []
        self.current_pris_iteration = 0

    def setup_pris_options(self, pris_options):
        """
        setup the pris options with input options. Otherwise, it will be the default options

        :param pris_options: the input pris_options.
        :return:

        """
        self.refinement_number = pris_options.refinement_number  # total number of refinement
        self.refinement_factor_list = pris_options.refinement_factor_list  # list of refinement factors
        self.lasso_solver = pris_options.lasso_solver  # choice of lasso solver. Only available one currently.
        self.dual_channel = pris_options.dual_channel

    def setup_lbreg(self, A, b, lbreg_opts):
        self.A = A
        self.b = b
        self.x = np.ndarray()
        self.lbreg_opts = lbreg_opts

    def prepare_next_iteration(self, pris_options, psf_model, coordinates_refiner):
        # based on the iteration number, do the following:
        # update location coordinates candidate voxels.
        coordinates_refiner(self, pris_options)
        # update A using all the candidate voxels, and PSF model.
        self.A = psf_model(self)

    def pris_go(self, lasso_solver):
        self.current_pris_iteration = self.current_pris_iteration + 1  # Pris iteration number tag increase by one.
        self.x = lasso_solver(self.A, self.b)


class PrisParaOpts:
    """
    I want to force objects of this class to have fixed number of attributes.

    """
    def __init__(self):
        self.type = 'Parameters'
        self.refinement_number = 5  # default total number of refinement
        self.refinement_factor_list = [1, 2, 2, 2, 2]  # default list of refinement factors
        self.lasso_solver = 'lbreg'  # default choice of lasso solver. Only available one currently.
        self.dual_channel = 'Yes'  # default setting of the number of observation channels.


class Lbreg:
    """ linearized bregman iteration.

            Keyword arguments:
            A -- the sensing matrix. data type: rank-2 tensor.
            f -- the input observation vector
            lbreg_opts -- options and parameters for the linearized bregman iteration.
    """
    def __init__(self, A, b, lbreg_opts):
        self.type = 'LASSO Solver'
        self.kick_switch = lbreg_opts.kick_switch
        self.A = A
        self.b = b
        self.v = []
        self.u = []
        "self.u = zeros(size(A(:,1))); " " -- abuse of MATLAB syntax."
    def update_v(self):
        """
          seems like the A'(Au - f) can be calculated in the form of the backpropagation of (Au - f)in PyTorch.
          the forward model will just be convolution, with pooling.
          :return: updated v
        """
    def update_u(self):
        """

        :return: updated u
        """

class LbregParaOpts:
    def __init__(self):
        self.type = 'Options'