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
import yaml

class RiverModule(object):

    def __init__(self):
        self._values = {}
        self._shape = (1000, 500)
        self._spacing = (10, 10)
        self._n0 = 100
        self._max_rand = 0.00001
        self._nslope = 0.0001
        self._spinup = 0
        self._dt_day = 73
        self._Initial_SL = 0
        self._SLRR_m = 0.015
        self._IRR_m = 0.005
        self._ch_width = 2000.
        self._ch_depth = 5.0
        self._init_cut_frac = 1
        self._nu = 10000
        self._super_ratio = 1
        self._short_path = 1
        self._time_max = 1
        self._time = 0.
        self._dt = 0.
        self._riv_x = None
        self._riv_y = None
        self._sed_flux = 0.
        self._shoreline = None

    @property
    def time(self):
        """Current model time."""
        return self._time

    @property
    def time_step(self):
        """Model time step."""
        return self._dt

    @time_step.setter
    def time_step(self, time_step):
        """Set model time step."""
        self._dt = time_step

    @property
    def river_x_coordinates(self):
        return self._riv_x

    @property 
    def river_y_coordinates(self):
        return self._riv_y

    @property 
    def sediment_flux(self):
        return self._sed_flux

    @property
    def avulsions(self):
        return self._avulsions

    @shoreline.setter
    def shoreline(self, shoreline):
        self._shoreline = shoreline

    @property
    def savefiles(self):
        return self._savefiles

    @classmethod
    def params_from_file(self, fname):
        """ create a RiverModule object from a file-like object. """

        with open(fname, 'r') as fp:
            params = yaml.load(fp)

        # Spatial parameters
        length, width = params['shape']
        dx_km, dy_km = params['spacing']
        self._L = length * 1000       # convert to meters
        self._W = width * 1000
        self._dx = dx_km * 1000
        self._dy = dy_km * 1000
        self._n0 = params['n0']
        self._max_rand = params['max_rand']
        self._nslope = params['nslope']

        self._imax = self._L/self._dx + 1
        self._jmax = self._W/self._dy + 1
        self._x = np.zeros((self._imax, self._jmax))   # longitudinal space
        self._y = np.zeros((self._imax, self._jmax))   # transverse space
        self._n = np.zeros((self._imax, self._jmax))   # eta, elevation
        self._dn_rc = np.zeros((self._imax))       # change in elevation along river course
        self._dn_fp = np.zeros((self._imax, self._jmax))     # change in elevation due to floodplain dep
        self._riv_x = [0]             # defines first x river locations
        self._riv_y = [self._W/2]          # defines first y river locations
        self._profile = np.zeros((self._imax))  # elevation profile of river course
        self._avulsions = [(0, 0, 0, 0, 0, 0)]    # initializes timestep/avulsions array

        # Time parameters
        self._dt = (params['dt_day'] *60*60*24)     # convert timestep to seconds
        self._time_max_s = (params['time_max'] * 31536000)  # length of model run in seconds
        self._spinup_s = (params['spinup'] * 31536000)  # length of spinup in seconds
        self._kmax = self._spinup_s/self._dt + self._time_max_s/self._dt + 1  # max number of timesteps
        self._save_after = self.spinup_s/self._dt        # save files after this point
        self._time = 0.
        self._k = 0

        # Sea level and subsidence parameters
        self._SL = [params['Initial_SL']]                   # initializes SL array
        self._SLRR = (params['SLRR_m'] / 31536000) * self._dt  # sea level rise rate in m/s per timestep
        self._IRR = (params['IRR_m'] / 31536000) * self._dt    # inlet rise rate in m/s per timestep
        self._shoreline = None  

        # River parameters
        self._nu = params['nu']
        self._init_cut = params['init_cut_frac'] * params['ch_depth']
        self._super_ratio = params['super_ratio']
        self._short_path = params['short_path']

        # Floodplain and wetland characteristics
        self._WL_Z = params['WL_Z']
        self._WL_dist = params['WL_dist']
        self._blanket_rate = (params['blanket_rate_m'] / 31536000) * self._dt    # blanket deposition in m/s
        self._splay_dep = (params['splay_dep_m'] / 31536000) * self._dt       # splay deposition in m/s
        self._splay_type = params['splay_type']

        # Saving information
        self._savefiles = params['savefiles']
        self._savespacing = params['savespacing']

        # Initialize elevation grid
        for i in range(self._imax):
            for j in range(self._jmax):
                self._x[i][j] = i * self._dx
                self._y[i][j] = j * self._dy
                self._n[i][j] = self._n0 - (self._nslope * float(self._x[i][j]) \
                          + self._max_rand * np.random.rand())
                j += 1
            i += 1

        # Determine initial river course
        self._riv_x, self._riv_y = steep_desc.find_course(self._dx, self._dy,
                                                          self._imax, self._jmax, self._n,
                                                          self._riv_x, self._riv_y)

        # downcut into new river course by amount determined by init_cut
        self._n = downcut.cut_init(self._dx, self._dy, self._riv_x, self._riv_y, self._n,
                                  self._init_cut)

        # smooth initial river course elevations using linear diffusion equation
        self._n, self._dn_rc = diffuse.smooth_rc(self._dx, self._dy, self._nu, self._dt,
                                                 self._riv_x, self._riv_y, self._n,
                                                 self._nslope)

        # Determine initial river profile
        self._profile = prof.make_profile(self._dx, self._dy, self._n, self._riv_x,
                                         self._riv_y, self._profile)

        # make directories and save initial condition files
        if self.savefiles == 1:
            # os.mkdir("run" + str(run_num) + "_out")
            os.mkdir("elev_grid")
            os.mkdir("riv_course")
            os.mkdir("profile")
            os.mkdir("dn_fp")
        #   saves initial conditions
        #    np.savetxt('elev_grid/elev_0.out', n, fmt='%f')
        #    np.savetxt('riv_course/riv_0.out', zip(riv_x, riv_y), fmt='%i')
        #    np.savetxt('profile/prof_0.out', profile, fmt='%f')

        # ### need to add self.var_units (talk to brad) ###
        return params

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
        self._riv_x, self._riv_y, self._loc, self._SEL, self._SER, self._n, \
            self._dn_fp, self._avulsion_type, self._length_new_sum, self._length_old \
            = avulse.find_avulsion(self._dx, self._dy, self._imax, self._jmax,
                                    self._riv_x, self._riv_y, self._n, self._super_ratio,
                                    self._current_SL, self._ch_depth, self._short_path,
                                    self._dn_fp, self._splay_type, self._splay_dep)

        # save timestep and avulsion location if there was one
        if len(self._loc) != 0:
            self._avulsions = self._avulsions + [(self._k*(self._dt/86400),
                        self._loc[-1], self._avulsion_type, self._length_old,
                        self._length_new_sum, self._current_SL)]
        
        # raise first two rows by inlet rise rate (subsidence)
        self._n[0][:] = self._n[0][:] + (self._IRR)
        self._n[1][:] = self._n[1][:] + (self._IRR)

        # change elevations according to sea level rise (SLRR)
        ### needs to be changed to subtracting elevation once coupled ###
        self._n, self._rc_flag = SLR.elev_change(self._imax, self._jmax,
                                                 self._current_SL, self._n, self._riv_x,
                                                 self._riv_y, self._ch_depth, self._dx,
                                                 self._dy)

        # smooth river course elevations using linear diffusion equation
        self._n, self._dn_rc = diffuse.smooth_rc(self._dx, self._dy, self._nu,
                                                 self._dt, self._riv_x, self._riv_y,
                                                 self._n, self._nslope)

        # Floodplain sedimentation
        self._n, self._dn_fp = FP.dep_blanket(self._dy, self._dx, self._imax,
                                              self._jmax, self._current_SL,
                                              self._blanket_rate, self._n, self._riv_x,
                                              self._riv_y, self._ch_depth)

        # Wetland sedimentation
        ### no wetlands in first version of coupling to CEM ###
        self._n, self._dn_fp = FP.wetlands(self._dx, self._dy, self._imax, self._jmax,
                                           self._current_SL, self._WL_Z, self._WL_dist,
                                           self._n, self._riv_x, self._riv_y, self._x,
                                           self._y, self._dn_fp)

        # calculate sediment flux
        self._sed_flux = flux.calc_qs(self._nu, self._riv_x, self._riv_y, self._n,
                                     self._dx, self._dy, self._dt)

        # create a river profile array
        self._profile = prof.make_profile(self._dx, self._dy, self._n, self._riv_x,
                                          self._riv_y, self._profile)

        #self._riv_mouth = [self._riv_x[-1], self._riv_y[-1]]

        self._k += 1
        self._time += self._dt

        # save files
        if model._savefiles == 1:
            if model._k >= model._save_after:
                if model._k % model._savespacing == 0:
                    np.savetxt('elev_grid/elev_' + str(model._k*model._dt/86400 
                                - model._save_after) + '.out', model._n, fmt='%.6f')
                    np.savetxt('riv_course/riv_' + str(model._k*model._dt/86400
                                - model._save_after) + '.out',
                                zip(model._riv_x, model._riv_y), fmt='%i')
                    np.savetxt('profile/prof_' + str(model._k*model._dt/86400
                                - model._save_after) + '.out',
                                model._profile, fmt='%.6f')
                    np.savetxt('dn_fp/dn_fp_' + str(model._k*model._dt/86400
                                - model._save_after) + '.out',
                                model._dn_fp, fmt='%.6f')

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

def main ():
    model = RiverModule()
    model.initialize('input.yaml')

    while model._k < model._kmax:
        model.advance_in_time()

        # # save files
        # if model._savefiles == 1:
        #     if model._k >= model._save_after:
        #         if model._k % model._savespacing == 0:
        #             np.savetxt('elev_grid/elev_' + str(model._k*model._dt/86400 
        #                         - model._save_after) + '.out', model._n, fmt='%.6f')
        #             np.savetxt('riv_course/riv_' + str(model._k*model._dt/86400
        #                         - model._save_after) + '.out',
        #                         zip(model._riv_x, model._riv_y), fmt='%i')
        #             np.savetxt('profile/prof_' + str(model._k*model._dt/86400
        #                         - model._save_after) + '.out',
        #                         model._profile, fmt='%.6f')
        #             np.savetxt('dn_fp/dn_fp_' + str(model._k*model._dt/86400
        #                         - model._save_after) + '.out',
        #                         model._dn_fp, fmt='%.6f')

    model.finalize()

if __name__ == '__main__':
    main()


