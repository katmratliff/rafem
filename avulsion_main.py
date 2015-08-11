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

class RiverModule(object):

    def __init__(self):
#        self._shape = (1000, 500)
#        self._spacing = (10, 10)
#        self._n0 = 100
#        self._max_rand = 0.00001
#        self._nslope = 0.0001
#        self._spinup = 0
#        self._dt_day = 73
#        self._Initial_SL = 0
#        self._SLRR_m = 0.015
#        self._IRR_m = 0.005
#        self._ch_width = 2000.
#        self._ch_depth = 5.0
#        self._init_cut_frac = 1
#        self._nu = 10000
#        self._super_ratio = 1
#        self._short_path = 1
#        self._time_max = 1
#        self._time = 0.
#        self._dt = 0.
#        self._riv_x = 0
#        self._riv_y = 0
#        self._sed_flux = 0.
        #self._shoreline = None
        pass

    @property
    def time(self):
        """Current model time (converted from seconds to days)."""
        return (self._time/86400)

    @property
    def time_step(self):
        """Model time step (converted from seconds to days)."""
        return (self._dt/86400)

    @time_step.setter
    def time_step(self, time_step):
        """Set model time step (time_step is in days)."""
        self._dt = (time_step*86400)

    @property
    def river_x_coordinates(self):
        return self._riv_i * self._dx

    @property 
    def river_y_coordinates(self):
        return self._riv_j * self._dy

    @property 
    def elevation(self):
        return self._n

    @property 
    def sediment_flux(self):
        return [self._sed_flux]

    #@property
    #def avulsions(self):
    #    return self._avulsions

#    @shoreline.setter
#    def shoreline(self, shoreline):
#        self._shoreline = shoreline

    @property
    def savefiles(self):
        return self._savefiles

    @classmethod
    def from_path(cls, fname):
        """Create a RiverModule object from a file-like object."""

        avulsion = cls()
        avulsion._init_from_dict(read_params_from_file(fname))

        return avulsion

    def _init_from_dict(self, params):
        # Spatial parameters
        length, width = params['shape']
        dx_km, dy_km = params['spacing']
        self._L = length * 1000.       # convert to meters
        self._W = width * 1000.
        self._dx = dx_km * 1000.
        self._dy = dy_km * 1000.
        self._n0 = params['n0']
        self._max_rand = params['max_rand']
        self._nslope = params['nslope']

        n_rows = int(self._L // self._dx + 1)
        n_cols = int(self._W // self._dy + 1)

        # Initialize elevation grid
        # transverse and longitudinal space
        self._y, self._x = np.meshgrid(np.arange(n_cols) * self._dy,
                                      np.arange(n_rows) * self._dx)
        # eta, elevation
        self._n = self._n0 - (self._nslope * self._x +
                              np.random.rand(n_rows, n_cols) * self._max_rand)

        #self._dn_rc = np.zeros((self._imax))       # change in elevation along river course
        self._dn_fp = np.zeros_like(self._n)     # change in elevation due to floodplain dep
        self._riv_i = np.zeros(1, dtype=np.int) # defines first x river locations
        self._riv_j = np.zeros(1, dtype=np.int) # defines first y river locations
        self._riv_j[0] = int(self._W / self._dx * .5)

        # Time parameters
        self._dt = params['dt_day'] * 60 * 60 * 24     # convert timestep to seconds
        self._time = 0.
        self._k = 0

        # Sea level and subsidence parameters
        self._SL = [params['Initial_SL']]                   # initializes SL array
        self._SLRR = (params['SLRR_m'] / 31536000) * self._dt  # sea level rise rate in m/s per timestep
        self._IRR = (params['IRR_m'] / 31536000) * self._dt    # inlet rise rate in m/s per timestep

        # River parameters
        self._nu = params['nu']
        self._init_cut = params['init_cut_frac'] * params['ch_depth']
        self._super_ratio = params['super_ratio']
        self._short_path = params['short_path']
        self._ch_width = params['ch_width']
        self._ch_depth = params['ch_depth']

        # Floodplain and wetland characteristics
        self._WL_Z = params['WL_Z']
        self._WL_dist = params['WL_dist']
        self._blanket_rate = (params['blanket_rate_m'] / 31536000) * self._dt    # blanket deposition in m/s
        self._splay_dep = (params['splay_dep_m'] / 31536000) * self._dt       # splay deposition in m/s
        self._splay_type = params['splay_type']

        self._sed_flux = 0.

        # Saving information
        #self._savefiles = params['savefiles']
        #self._savespacing = params['savespacing']

        self._riv_i, self._riv_j = steep_desc.find_course(self._n, self._riv_i,
                                                          self._riv_j)

        # downcut into new river course by amount determined by init_cut
        downcut.cut_init(self._riv_i, self._riv_j, self._n, self._init_cut)

        # smooth initial river course elevations using linear diffusion equation
        diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt,
                          self._riv_i, self._riv_j, self._n)

    def advance_in_time(self):
        """ Update avulsion model one time step. """

        # begin time loop and main program
        # for k in range(kmax):

        # determine current sea level
        self._SL = self._SL + [self._k * self._SLRR]
        self._current_SL = self._SL[-1]

        ### future work: SLRR can be a vector to change rates ###

        # determine if there is an avulsion & find new path if so
        ### need to change this to look for shoreline after coupling ###
        ### (instead of looking for sea level)
        self._riv_i, self._riv_j = avulse.find_avulsion(
             self._riv_i, self._riv_j, self._n,
             self._super_ratio, self._current_SL, self._ch_depth,
             self._short_path, self._splay_type, self._splay_dep)

        #assert(self._riv_i[-1] != 0)

        # save timestep and avulsion location if there was one
        #if len(loc) != 0:
        #    self._avulsions = self._avulsions + [(self._k*(self._dt/86400),
        #                loc[-1], avulsion_type, length_old,
        #                length_new_sum, self._current_SL)]
        
        # raise first two rows by inlet rise rate (subsidence)
        self._n[:2, :] += self._IRR

        # change elevations according to sea level rise (SLRR)
        ### needs to be changed to subtracting elevation once coupled ###
        SLR.elev_change(self._current_SL, self._n, self._riv_i,
                        self._riv_j, self._ch_depth)

        # smooth river course elevations using linear diffusion equation
        diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt,
                          self._riv_i, self._riv_j, self._n)

        # Floodplain sedimentation
        FP.dep_blanket(self._current_SL, self._blanket_rate, self._n,
                       self._riv_i, self._riv_j, self._ch_depth)

        # Wetland sedimentation
        ### no wetlands in first version of coupling to CEM ###
        FP.wetlands(self._current_SL, self._WL_Z, self._WL_dist * self._dy,
                    self._n, self._riv_i, self._riv_j, self._x, self._y)

        # calculate sediment flux
        self._sed_flux = flux.calc_qs(self._nu, self._riv_i,
                                      self._riv_j, self._n,
                                      self._dx, self._dy, self._dt)

        # create a river profile array
        #self._profile = prof.make_profile(self._dx, self._dy, self._n,
        #                                  np.array(self._riv_x),
        #                                  np.array(self._riv_y),
        #                                  self._profile)

        #self._riv_mouth = [self._riv_x[-1], self._riv_y[-1]]

        self._k += 1
        self._time += self._dt
