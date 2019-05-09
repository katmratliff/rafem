#! /usr/local/bin/python

import os
import numpy as np
import errno

import yaml

from . import (
    steep_desc,
    avulse,
    diffuse,
    prof,
    SLR,
    FP,
    downcut,
    flux,
    subside,
    avulsion_utils,
)
from .avulsion_utils import read_params_from_file


_SECONDS_PER_YEAR = 31536000.
_SECONDS_PER_DAY = 86400.


def make_empty_file(path):
    """Create an empty file.

    Create an empty file along with all of its parent folders,
    if necessary. Note that if the file already exists, it
    will be clobbered.

    Parameters
    ----------
    path : str
        Path to the file to create.
    """
    dirname = os.path.dirname(path)
    try:
        os.makedirs(dirname)
    except OSError as err:
        if err.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise

    with open(path, "w"):
        pass


class RiverModule(object):

    def __init__(
        self,
        rand_seed=1945,
        spacing=(0.1, 0.1),
        shape=(100, 120),
        n0=5.0,
        nslope=0.001,
        max_rand=0.1,
        dt_day=0.05,
        Initial_SL=0.0,
        SLRR_m=0.0,
        SubRate_m=0.0,
        Sub_Start=0,
        ch_width=10.0,
        ch_depth=1.0,
        ch_discharge=10.0,
        A=1.,
        c_f=0.01,
        C_0=1.,
        sed_sg=2.65,
        init_cut_frac=1.0,
        super_ratio=1.0,
        short_path=1,
        WL_Z=0.0,
        WL_dist=0,
        blanket_rate_m=0.0,
        fine_dep_frac=0.0,
        splay_type=2,
        saveavulsions=False,
        savecourseupdates=False,
    ):
        """The RAFEM.

        Parameters
        ----------
        rand_seed : int
            Seed for random number generator.
        spacing : tuple of float
            Row (dy) and column (dx) spacing [m].
        shape : tuple of float
            Size of grid as *(n_rows, n_columns)*.
        n0 : float
            Upstream elevation [m].
        nslope : float
            Initial landscape slope [-].
        max_rand : float
            Multiply by slope for max height of a random perturbation.
        dt_day : float
            Time step [days].
        Initial_SL : float
            Initial sea level [m].
        SLRR_m : float
            Sea level rise rate [m / yr].
        SubRate_m: float
            Subsidence rate [m / yr].
        Sub_Start : int
            Row where subsidence starts
        ch_width : float
            Characteristic channel width [m].
        ch_depth : float
            Characteristic channel depth [m].
        ch_discharge : float
            Long-term averaged discharge [m^3 / s].
        A : float
            River-dependent const (1 for meandering; 1.4 for braided).
        c_f : float
            Drag coefficient.
        C_0 : float
            Sediment concentration on bed.
        sed_sg : float
            Sediment specific gravity
        init_cut_frac : float
            Initial cut of the channel into land surface.
        super_ratio : float
            Normalized SE ratio to trigger avulsion.
        short_path : int
            Flag for using shortest path to complete avulsion.
        WL_Z : float
            Elevation that wetlands maintain above SL [m].
        WL_dist : float
            Cell distance beyond channel that wetlands exist.
        blanket_rate_m : float
            "Blanket" deposition rate [m].
        fine_dep_frac : float
            Fraction of channel deposit for adjacent fine deposition.
        splay_type : {0, 1, 2}
            Size/type of splay. Values of *splay_type* correspond to:
            *  0: no splay deposition
            *  1: just the first failed avulsion river cell
            *  2: first failed cell and adjacent cells
        saveavulsions : bool
            Flag for saving avulsion info.
        savecourseupdates : bool
            Flag for saving course updates
        """
        self._params = locals().copy()
        self._params.pop("self")

        np.random.seed(rand_seed)

        # Spatial parameters
        self._dy = spacing[0] * 1000.
        self._dx = spacing[1] * 1000.

        n_rows = int(shape[0])
        n_cols = int(shape[1])

        # Initialize elevation grid
        # transverse and longitudinal space
        self._x, self._y = np.meshgrid(np.arange(n_cols) * self._dx,
                                       np.arange(n_rows) * self._dy)
        # eta, elevation
        self._slope = nslope
        self._max_rand = max_rand * nslope
        self._n = n0 - (self._slope * self._y +
                        np.random.rand(n_rows, n_cols) * self._max_rand)
        self._n -= 0.05

        # self._dn_rc = np.zeros((self._imax))       # change in elevation along river course
        # self._dn_fp = np.zeros_like(self._n)     # change in elevation due to floodplain dep

        self._riv_i = np.zeros(1, dtype=np.int) # defines first x river locations
        self._riv_j = np.zeros(1, dtype=np.int) # defines first y river locations
        self._riv_j[0] = self._n.shape[1] / 2

        # Time parameters
        self._dt = dt_day * _SECONDS_PER_DAY # convert timestep to seconds
        self._time = 0.

        # Sea level and subsidence parameters
        self._SL = Initial_SL  # starting sea level
        self._SLRR = SLRR_m / _SECONDS_PER_YEAR * self._dt  # sea level rise rate in m (per timestep)
        self._SubRate = SubRate_m / _SECONDS_PER_YEAR * self._dt  # subsidence rate in m (per timestep)
        self._SubStart = Sub_Start  # row where subsidence begins

        # River parameters
        self._nu = ((8. * (ch_discharge / ch_width) * A
                    * np.sqrt(c_f)) / (C_0 * (sed_sg - 1)))
        ### NEED TO REDO DIFFUSE.PY TO HAVE SIGN OF NU CORRECT (NEG) ABOVE ###
        init_cut = init_cut_frac * ch_depth
        self._super_ratio = super_ratio
        self._short_path = short_path
        self._ch_depth = ch_depth

        # Floodplain and wetland characteristics
        self._WL_Z = WL_Z
        self._WL_dist = WL_dist
        self._blanket_rate = (blanket_rate_m / _SECONDS_PER_YEAR) * self._dt    # blanket deposition in m
        self._splay_type = splay_type
        self._frac_fines = fine_dep_frac

        self._sed_flux = 0.
        self._splay_deposit = np.zeros_like(self._n)

        # Saving information
        self._saveavulsions = saveavulsions
        self._saveupdates = savecourseupdates
        self._save_splay_deposit = False

        if self._saveupdates:
            make_empty_file(self._saveupdates)
        if self._saveavulsions:
            make_empty_file(self._saveavulsions)
        if self._save_splay_deposit:
            make_empty_file(self._savesplay_deposit)

        self._riv_i, self._riv_j = steep_desc.find_course(self._n, self._riv_i, self._riv_j,
                                                          len(self._riv_i), self._ch_depth,
                                                          sea_level=self._SL)

        # downcut into new river course by amount determined by init_cut
        downcut.cut_init(self._riv_i, self._riv_j, self._n, init_cut)

        # smooth initial river course elevations using linear diffusion equation
        diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt, self._ch_depth,
                          self._riv_i, self._riv_j, self._n, self._SL, self._slope)

        # initial profile
        self._profile = self._n[self._riv_i, self._riv_j]

        self._profile = self._n[self._riv_i, self._riv_j]

    def to_yaml(self):
        return yaml.dump(self._params)

    @property
    def time(self):
        """Current model time (converted from seconds to days)."""
        return self._time / _SECONDS_PER_DAY

    @property
    def time_step(self):
        """Model time step (converted from seconds to days)."""
        return self._dt / _SECONDS_PER_DAY

    @time_step.setter
    def time_step(self, time_step):
        """Set model time step (time_step is in days)."""
        self._dt = time_step * _SECONDS_PER_DAY

    @property
    def grid_shape(self):
        return self._n.shape

    @property
    def grid_spacing(self):
        return (self._dy, self._dx)

    @property
    def river_x_coordinates(self):
        return self._riv_j * self._dx

    @river_x_coordinates.setter
    def river_x_coordinates(self, new_coords):
        self._riv_j = np.asarray(new_coords) / self._dx

    @property 
    def river_y_coordinates(self):
        return self._riv_i * self._dy

    @river_y_coordinates.setter
    def river_y_coordinates(self, new_coords):
        self._riv_i = np.asarray(new_coords) / self._dy

    @property
    def sea_level(self):
        return self._SL

    @property 
    def elevation(self):
        return self._n
        # return self._n + self.sea_level

    @elevation.setter
    def elevation(self, new_elev):
        """Set the land surface elevation."""
        self._elevation[:] = new_elev
        # self._elevation[:] = new_elev - self.sea_level

    @property 
    def sediment_flux(self):
        return self._sed_flux

    @property
    def avulsions(self):
        return self._avulsion_info

    @property
    def profile(self):
        self._profile[-1] = self._SL - self._ch_depth
        return self._profile

#    @shoreline.setter
#    def shoreline(self, shoreline):
#        self._shoreline = shoreline

    @classmethod
    def from_path(cls, fname):
        """Create a RiverModule object from a file-like object."""
        if fname:
            params = read_params_from_file(fname)
        else:
            params = dict()
        return cls(**params)

    def advance_in_time(self):
        """ Update avulsion model one time step. """
        # if (self._time / _SECONDS_PER_YEAR) > 2000:
        #     self._SLRR = 0.01 / _SECONDS_PER_YEAR * self._dt

        self._riv_i, self._riv_j, self._course_update = steep_desc.update_course(
            self._n, self._riv_i, self._riv_j, self._ch_depth, self._slope,
            sea_level=self._SL, dx=self._dx, dy=self._dy)

        self._n = avulsion_utils.fix_elevations(self._n, self._riv_i, self._riv_j,
            self._ch_depth, self._SL, self._slope, self._dx, self._max_rand, self._SLRR)

        """ Save every time the course changes? """
        if self._saveupdates and self._course_update > 0:
            with open(self._saveupdates,'a') as file:
                file.write("%.5f %i \n" % ((self._time / _SECONDS_PER_YEAR),
                    self._course_update))

        """ determine if there is an avulsion & find new path if so """
        (self._riv_i, self._riv_j), self._avulsion_type, self._loc, self._avulse_length, \
         self._path_diff, self._splay_deposit = avulse.find_avulsion(self._riv_i,
            self._riv_j, self._n, self._super_ratio, self._SL, self._ch_depth,
            self._short_path, self._splay_type, self._slope,
            self._splay_deposit, self._nu, self._dt, dx=self._dx, dy=self._dy)

        """ Save avulsion record. """
        if self._saveavulsions and self._avulsion_type > 0:
            with open(self._saveavulsions,'a') as file:
                file.write("%.5f %i %i %.5f %.5f\n" % ((self._time / _SECONDS_PER_YEAR),
                    self._avulsion_type, self._loc, self._avulse_length, self._path_diff))

        """ Save crevasse splay deposits. """        
        if self._save_splay_deposit and (self._splay_deposit.sum() > 0):
            np.savetxt(self._save_splay_deposit, self._splay_deposit, '%.8f')

        # need to fill old river channels if coupled to CEM
        if (self._avulsion_type == 1) or (self._avulsion_type == 2):
            self._n = avulsion_utils.fix_elevations(self._n, self._riv_i, self._riv_j,
                self._ch_depth, self._SL, self._slope, self._dx, self._max_rand, self._SLRR)

        #assert(self._riv_i[-1] != 0)

        """ change elevations according to sea level rise (SLRR)
        (if not coupled -- this occurs in coupling script otherwise) """
        # SLR.elev_change(self._SL, self._n, self._riv_i,
        #                 self._riv_j, self._ch_depth, self._SLRR)

        """ smooth river course elevations using linear diffusion equation """
        self._dn_rc = diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt, self._ch_depth,
                          self._riv_i, self._riv_j, self._n, self._SL, self._slope)

        """ Floodplain sedimentation (use one or the other) """
        #-------------------------------------------------------
        ### Deposit blanket across entire subaerial domain: ###
        # FP.dep_blanket(self._SL, self._blanket_rate, self._n,
        #                self._riv_i, self._riv_j, self._ch_depth)

        ### Deposit fines adjacent to river channel: ###
        FP.dep_fines(self._n, self._riv_i, self._riv_j, self._dn_rc, self._frac_fines,
                     self._SL)
        #-------------------------------------------------------

        """ Wetland sedimentation """
        ### no wetlands in first version of coupling to CEM ###
        # FP.wetlands(self._SL, self._WL_Z, self._WL_dist * self._dy,
        #             self._n, self._riv_i, self._riv_j, self._x, self._y)

        """ Subsidence """
        subside.linear_subsidence(self._n, self._riv_i, self._riv_j, self._ch_depth,
                                  self._SubRate, self._SubStart, self._SL)

        """ calculate sediment flux at the river mouth """
        self._sed_flux = flux.calc_qs(self._nu, self._riv_i, self._riv_j,
                                      self._n, self._SL, self._ch_depth,
                                      self._dx, self._dy, self._dt, self._slope)

        self._profile = self._n[self._riv_i, self._riv_j]

        # Update time
        self._time += self._dt
        # update sea level
        self._SL += self._SLRR
