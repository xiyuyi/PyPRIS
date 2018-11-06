import numpy as np

class LbregSolver:
    """ linearized bregman iteration.

            Keyword arguments:
            A -- the sensing matrix. data type: rank-2 tensor.
            f -- the input observation vector
            lbreg_opts -- options and parameters for the linearized bregman iteration.
    """
    def __init__(self, A, b, lbreg_opts):
        self.type = 'sparse recovery solver'
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

class LbregOpts:
    def __init__(self):
        self.type = 'options for linearized bregman iteration'

