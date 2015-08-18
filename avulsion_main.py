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
    def river_x_coordinates(self):
        return self._riv_j * self._dx

    @property 
    def river_y_coordinates(self):
        return self._riv_i * self._dy

    @property 
    def elevation(self):
        return self._n

    @property 
    def sediment_flux(self):
        return self._sed_flux

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
    def params_from_file(self, fname):
        """ create a RiverModule object from a file-like object. """

        params = read_params_from_file(fname)
        #params = yaml.load(fname)

        # Spatial parameters
        self._dy = params['spacing'][0] * 1000.
        self._dx = params['spacing'][1] * 1000.

        length, width = params['shape']
        n_rows = int(length * 1000 / self._dy + 1)
        n_cols = int(width * 1000 / self._dx + 1)

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
        self._dt = params['dt_day'] * 60. * 60. * 24. # convert timestep to seconds
        self._time = 0.

        # Sea level and subsidence parameters
        self._SL = params['Initial_SL'] # starting sea level
        self._SLRR = params['SLRR_m'] / _SECONDS_PER_YEAR * self._dt # sea level rise rate in m (per timestep)
        self._IRR = params['IRR_m'] / _SECONDS_PER_YEAR * self._dt # inlet rise rate in m

        # River parameters
        self._nu = params['nu'] / _SECONDS_PER_DAY
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

        # Saving information
        #self._savefiles = params['savefiles']
        #self._savespacing = params['savespacing']

        self._riv_i, self._riv_j = steep_desc.find_course(self._n, self._riv_i,
                                                          self._riv_j)

        # downcut into new river course by amount determined by init_cut
        downcut.cut_init(self._riv_i, self._riv_j, self._n, init_cut)

        # smooth initial river course elevations using linear diffusion equation
        diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt,
                          self._riv_i, self._riv_j, self._n)

    def advance_in_time(self):
        """ Update avulsion model one time step. """

        ### future work: SLRR can be a vector to change rates ###

        # determine if there is an avulsion & find new path if so
        ### need to change this to look for shoreline after coupling ###
        ### (instead of looking for sea level)
        self._riv_i, self._riv_j = avulse.find_avulsion(
             self._riv_i, self._riv_j, self._n,
             self._super_ratio, self._SL, self._ch_depth,
             self._short_path, self._splay_type, self._splay_dep, dx=self._dx,
             dy=self._dy)

        #assert(self._riv_i[-1] != 0)

        # save timestep and avulsion location if there was one
        #if len(loc) != 0:
        #    self._avulsions = self._avulsions + [(self._k*(self._dt/86400),
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

        # Update sea level
        self._SL += self._SLRR
        self._time += self._dt

        # save files
        if self._savefiles == 1:
            if self._k >= self._save_after:
                if self._k % self._savespacing == 0:
                    np.savetxt('elev_grid/elev_' + str(self._k*self._dt/86400 
                                - self._save_after) + '.out', self._n, fmt='%.6f')
                    np.savetxt('riv_course/riv_' + str(self._k*self._dt/86400
                                - self._save_after) + '.out',
                                zip(self._riv_x, self._riv_y), fmt='%i')
                    np.savetxt('profile/prof_' + str(self._k*self._dt/86400
                                - self._save_after) + '.out',
                                self._profile, fmt='%.6f')
                    np.savetxt('dn_fp/dn_fp_' + str(self._k*self._dt/86400
                                - self._save_after) + '.out',
                                self._dn_fp, fmt='%.6f')

        # print "sediment flux = %f" % self.sed_flux

    # def finalize(self):
    #     """Clean up & save avulsion file"""
        
    #     if self.savefiles == 1:
    #         np.savetxt('avulsions', self._avulsions, fmt='%i %i %i %.3f %.3f %.3f')
    #     pass

    # def get_value_ref(self, var_name):
    #     return self._values[var_name]

    # # def get_value(self, var_name):
    # #     return self.get_value_ref(var_name).copy()

    # def get_value_at_indices(self, var_name, indices):
    #     return self.get_value_ref(var_name).take(indices)

    # def set_value(self, var_name, src):
    #     val = self.get_value_ref(var_name)
    #     val[:] = src

    # def set_value_at_indices(self, var_name, src, indices):
    #     val = self.get_value_ref(var_name)
    #     val.flat[indices] = src

    # def get_component_name(self):
    #     """Name of the component."""
    #     return self._name

    # def get_input_var_names(self):
    #     """Get names of input variables."""
    #     return self._input_var_names

    # def get_output_var_names(self):
    #     """Get names of output variables."""
    #     return self._output_var_names

    # def get_current_time(self):
    #     return self._model.time/86400

    # def get_time_step(self):
    #     """Time step of model."""
    #     return self._model.dt/86400

# def main ():
#     model = RiverModule()
#     model.initialize('input.yaml')

#     while model._k < model._kmax:
#         model.advance_in_time()

#     model.finalize()

# if __name__ == '__main__':
#     main()


