#! /usr/local/bin/python
"""
Created on Wed Nov 12 09:28:51 2014

@author: kmratliff
"""
# from pylab import *
import os
import numpy as np
import steep_desc
import avulse
import diffuse
import prof
import SLR
import FP
import downcut
import flux
from avulsion_utils import read_params_from_file


_SECONDS_PER_YEAR = 31536000.
_SECONDS_PER_DAY = 86400.


class RiverModule(object):

    def __init__(self):
        pass

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

    @elevation.setter
    def elevation(self, new_elev):
        """Set the land surface elevation."""
        self._elevation[:] = new_elev

    @property 
    def sediment_flux(self):
        return self._sed_flux

    @property
    def avulsions(self):
        return self._avulsion_info

    @property
    def profile(self):
        return self._profile

#    @shoreline.setter
#    def shoreline(self, shoreline):
#        self._shoreline = shoreline

    @classmethod
    def from_path(cls, fname):
        """Create a RiverModule object from a file-like object."""

        avulsion = cls()
        avulsion._init_from_dict(read_params_from_file(fname))

        return avulsion

    def _init_from_dict(self, params):
        # Spatial parameters
        self._dy = params['spacing'][0] * 1000.
        self._dx = params['spacing'][1] * 1000.

        n_rows = int(params['shape'][0])
        n_cols = int(params['shape'][1])
        #n_rows = int(length * 1000 / self._dy + 1)
        #n_cols = int(width * 1000 / self._dx + 1)

        # Initialize elevation grid
        # transverse and longitudinal space
        self._x, self._y = np.meshgrid(np.arange(n_cols) * self._dx,
                                       np.arange(n_rows) * self._dy)
        # eta, elevation
        n0 = params['n0']
        max_rand = params['max_rand']
        slope = params['nslope']
        self._n = n0 - (slope * self._y +
                        np.random.rand(n_rows, n_cols) * max_rand)

        #self._dn_rc = np.zeros((self._imax))       # change in elevation along river course
        #self._dn_fp = np.zeros_like(self._n)     # change in elevation due to floodplain dep

        self._riv_i = np.zeros(1, dtype=np.int) # defines first x river locations
        self._riv_j = np.zeros(1, dtype=np.int) # defines first y river locations
        self._riv_j[0] = self._n.shape[1] / 2

        # Time parameters
        self._dt = params['dt_day'] * _SECONDS_PER_DAY # convert timestep to seconds
        self._time = 0.

        # Sea level and subsidence parameters
        self._SL = params['Initial_SL'] # starting sea level
        self._SLRR = params['SLRR_m'] / _SECONDS_PER_YEAR * self._dt # sea level rise rate in m (per timestep)
        self._IRR = params['IRR_m'] / _SECONDS_PER_YEAR * self._dt # inlet rise rate in m

        # River parameters
        self._nu = ((8. * (params['ch_discharge'] / params['ch_width']) * params['A']
                    * np.sqrt(params['c_f'])) / (params['C_0'] * (params['sed_sg'] - 1)))
        ### NEED TO REDO DIFFUSE.PY TO HAVE SIGN OF NU CORRECT (NEG) ABOVE ###
        #self._nu = params['nu'] / _SECONDS_PER_DAY
        init_cut = params['init_cut_frac'] * params['ch_depth']
        self._super_ratio = params['super_ratio']
        self._short_path = params['short_path']
        self._ch_depth = params['ch_depth']

        # Floodplain and wetland characteristics
        self._WL_Z = params['WL_Z']
        self._WL_dist = params['WL_dist']
        self._blanket_rate = (params['blanket_rate_m'] / _SECONDS_PER_YEAR) * self._dt    # blanket deposition in m
        self._splay_dep = (params['splay_dep_m'] / _SECONDS_PER_YEAR) * self._dt       # splay deposition in m
        self._splay_type = params['splay_type']

        self._sed_flux = 0.
        self._avulsion_info = np.zeros(3, dtype=np.float)

        # Saving information
        self._saveavulsions = params['saveavulsions']
        #self._savefiles = params['savefiles']
        #self._savespacing = params['savespacing']

        self._riv_i, self._riv_j = steep_desc.find_course(self._n, self._riv_i,
                                                          self._riv_j,
                                                          sea_level=self._SL)

        # downcut into new river course by amount determined by init_cut
        downcut.cut_init(self._riv_i, self._riv_j, self._n, init_cut)

        # smooth initial river course elevations using linear diffusion equation
        diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt,
                          self._riv_i, self._riv_j, self._n)

    def advance_in_time(self):
        """ Update avulsion model one time step. """

        ### future work: SLRR can be a vector to change rates ###

        old_len = len(self._riv_i)
        self._riv_i, self._riv_j = steep_desc.find_course(
            self._n, self._riv_i, self._riv_j, sea_level=self._SL)

        # determine if there is an avulsion & find new path if so
        ### need to change this to look for shoreline after coupling ###
        ### (instead of looking for sea level)
        (self._riv_i, self._riv_j), self._avulsion_type, self._loc = avulse.find_avulsion(
             self._riv_i, self._riv_j, self._n,
             self._super_ratio, self._SL, self._ch_depth,
             self._short_path, self._splay_type, self._splay_dep, dx=self._dx,
             dy=self._dy)

        if self._saveavulsions & self._avulsion_type > 0:
            new_info = (self._avulsion_type, self._time / _SECONDS_PER_YEAR, self._loc)
            self._avulsion_info = np.vstack([self._avulsion_info, new_info])

        #assert(self._riv_i[-1] != 0)

        # save timestep and avulsion location if there was one
        #if len(loc) != 0:
        #    self._avulsions = self._avulsions + [(self._time/_SECONDS_PER_DAY),
        #                loc[-1], avulsion_type, length_old,
        #                length_new_sum, self._SL)]
        
        # raise first two rows by inlet rise rate (subsidence)
        self._n[:2, :] += self._IRR

        # change elevations according to sea level rise (SLRR)
        ### needs to be changed to subtracting elevation once coupled ###
        SLR.elev_change(self._SL, self._n, self._riv_i,
                        self._riv_j, self._ch_depth)

        # smooth river course elevations using linear diffusion equation
        diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt,
                          self._riv_i, self._riv_j, self._n)

        # Floodplain sedimentation
        FP.dep_blanket(self._SL, self._blanket_rate, self._n,
                       self._riv_i, self._riv_j, self._ch_depth)

        # Wetland sedimentation
        ### no wetlands in first version of coupling to CEM ###
        FP.wetlands(self._SL, self._WL_Z, self._WL_dist * self._dy,
                    self._n, self._riv_i, self._riv_j, self._x, self._y)

        # calculate sediment flux
        self._sed_flux = flux.calc_qs(self._nu, self._riv_i,
                                      self._riv_j, self._n,
                                      self._dx, self._dy, self._dt)

        self._profile = self._n[self._riv_i, self._riv_j]

        # Update sea level
        self._SL += self._SLRR
        self._time += self._dt
