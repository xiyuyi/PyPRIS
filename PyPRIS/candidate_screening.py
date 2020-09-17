import numpy as np
from PyPRIS import *

class ProjFilter2D_Opts:
    def __init__(self, thres = 0.2):
        self.relative_lower_bound_2Dproj = thres

class RelativeValueFilter_Opts:
    def __init__(self):
        self.relative_lower_bound_FullPool = 0.2

class PercentageVoxelCounts_Opts:
    def __init__(self):
        self.percentage_voxels_to_keep = 0.8
