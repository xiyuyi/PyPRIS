import numpy as np
import torch


class PRIS:
    """

    """
    def __init__(self):
        self.type = 'PRIS solver'
        self.A = []  # sensing matrix
        self.b = []  # observation
        self.x = torch.empty(1, 1, dtype=torch.float)
        self.sparse_recovery_opts = []
        self.current_pris_iteration = 0

    def pris_loader(self, inputs):

    def pris_refine(self, pris_options, psf_model, selector, refiner):
        # based on the iteration number, do the following:
        # identify the selected voxels
        selector(self, pris_options)
        # refine the selected candidate voxels.
        refiner(self, pris_options)
        # update A using all the candidate voxels, and PSF model.
        self.A = psf_model(self)

    def pris_go(self, subsolver):
        # based on the refined sensing matrix, perform sparse-recovery
        self.current_pris_iteration = self.current_pris_iteration + 1  # Pris iteration number tag increase by one.
        self.x = subsolver(self.A, self.b)


class PrisOpts:
    """
    I want to force objects of this class to have fixed number of attributes.

    """
    def __init__(self):
        self.type = 'Parameters'
        self.refinement_number = 5  # default total number of refinement
        self.refinement_factor_list = [1, 2, 2, 2, 2]  # default list of refinement factors
        self.subsolver = 'lbreg'  # default choice of lasso solver. Only available one currently.
        self.dual_channel = 'Yes'  # default setting of the number of observation channels.

class Selector:
    """

    """
    def __init__(self, pris_options, ):
        self.type = 'PRIS selector'

