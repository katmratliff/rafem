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
from bmi import Bmi

class River_Module(Bmi):
    _name = 'Avulsion module'
    _input_var_names = ('sea_shoreline')
    # not sure what's the most appropriate CSDMS river mouth location name?
    _output_var_names = ('river_mouth_location',
                         'channel_water_sediment~bedload__volume_flow_rate')

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
        self.nu = 10000
        self.super_ratio = 1
        self.short_path = 1
        self.time_max = 1

    def initialize(self, fname):
        """ initialize the avulsion module """

        params = read_params_from_file('input.yaml')

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
        self.n = np.zeros((self._imax, self._jmax))   # eta, elevation
        self.dn_rc = np.zeros((self._imax))       # change in elevation along river course
        self.dn_fp = np.zeros((self._imax, self._jmax))     # change in elevation due to floodplain dep
        self.riv_x = [0]             # defines first x river locations
        self.riv_y = [self._W/2]          # defines first y river locations
        self.profile = np.zeros((self._imax))  # elevation profile of river course
        self.avulsions = [(0, 0, 0, 0, 0, 0)]    # initializes timestep/avulsions array

        # Time parameters
        self.dt = (params['dt_day'] *60*60*24)     # convert timestep to seconds
        self.time_max_s = (params['time_max'] * 31536000)  # length of model run in seconds
        self.spinup_s = (params['spinup'] * 31536000)  # length of spinup in seconds
        self.kmax = self.spinup_s/self.dt + self.time_max_s/self.dt + 1  # max number of timesteps
        self.save_after = self.spinup_s/self.dt        # save files after this point
        self.time = 0.
        self.k = 0

        # Sea level and subsidence parameters
        self.SL = [params['Initial_SL']]                   # initializes SL array
        self.SLRR = (params['SLRR_m'] / 31536000) * self.dt  # sea level rise rate in m/s per timestep
        self.IRR = (params['IRR_m'] / 31536000) * self.dt    # inlet rise rate in m/s per timestep

        # River parameters
        self.nu = params['nu']
        self.init_cut = params['init_cut_frac'] * params['ch_depth']
        self.super_ratio = params['super_ratio']
        self.short_path = params['short_path']
        self.riv_mouth = None
        self.sed_flux = 0

        # Floodplain and wetland characteristics
        self.WL_Z = params['WL_Z']
        self.WL_dist = params['WL_dist']
        self.blanket_rate = (params['blanket_rate_m'] / 31536000) * self.dt    # blanket deposition in m/s
        self.splay_dep = (params['splay_dep_m'] / 31536000) * self.dt       # splay deposition in m/s
        self.splay_type = params['splay_type']

        # Saving information
        self.savefiles = params['savefiles']
        self.savespacing = params['savespacing']

        # Initialize elevation grid
        for i in range(self._imax):
            for j in range(self._jmax):
                self._x[i][j] = i * self._dx
                self._y[i][j] = j * self._dy
                self.n[i][j] = self._n0 - (self._nslope * float(self._x[i][j]) \
                          + self._max_rand * np.random.rand())
                j += 1
            i += 1

        # Determine initial river course
        self.riv_x, self.riv_y = steep_desc.find_course(self._dx, self._dy,
                                                        self._imax, self._jmax, self.n,
                                                        self.riv_x, self.riv_y)

        # downcut into new river course by amount determined by init_cut
        self.n = downcut.cut_init(self._dx, self._dy, self.riv_x, self.riv_y, self.n,
                                  self.init_cut)

        # smooth initial river course elevations using linear diffusion equation
        self.n, self.dn_rc = diffuse.smooth_rc(self._dx, self._dy, self.nu, self.dt,
                                               self.riv_x, self.riv_y, self.n, self._nslope)

        # Determine initial river profile
        self.profile = prof.make_profile(self._dx, self._dy, self.n, self.riv_x,
                                         self.riv_y, self.profile)

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

        self._values = {
            'river_mouth_location': self.riv_mouth,
            'channel_water_sediment~bedload__volume_flow_rate': self.sed_flux
        }

        ### need to add self.var_units (talk to brad) ###

    def update(self):
        """ Update avulsion model one time step. """

        # begin time loop and main program
        # for k in range(kmax):

        # determine current sea level
        self.SL = self.SL + [self.k * self.SLRR]
        self.current_SL = self.SL[-1]

        ### future work: SLRR can be a vector to change rates ###

        # determine if there is an avulsion & find new path if so
        ### need to change this to look for shoreline after coupling ###
        ### (instead of looking for sea level)
        self.riv_x, self.riv_y, self.loc, self.SEL, self.SER, self.n, \
            self.dn_fp, self.avulsion_type, self.length_new_sum, self.length_old \
            = avulse.find_avulsion(self._dx, self._dy, self._imax, self._jmax,
                                    self.riv_x, self.riv_y, self.n, self.super_ratio,
                                    self.current_SL, self._ch_depth, self.short_path,
                                    self.dn_fp, self.splay_type, self.splay_dep)

        # save timestep and avulsion location if there was one
        if len(self.loc) != 0:
            self.avulsions = self.avulsions + [(self.k*(self.dt/86400),
                        self.loc[-1], self.avulsion_type, self.length_old,
                        self.length_new_sum, self.current_SL)]
        
        # raise first two rows by inlet rise rate (subsidence)
        self.n[0][:] = self.n[0][:] + (self.IRR)
        self.n[1][:] = self.n[1][:] + (self.IRR)

        # change elevations according to sea level rise (SLRR)
        ### needs to be changed to subtracting elevation once coupled ###
        self.n, self.rc_flag = SLR.elev_change(self._imax, self._jmax,
                                               self.current_SL, self.n, self.riv_x,
                                               self.riv_y, self._ch_depth, self._dx,
                                               self._dy)

        # smooth river course elevations using linear diffusion equation
        self.n, self.dn_rc = diffuse.smooth_rc(self._dx, self._dy, self.nu,
                                               self.dt, self.riv_x, self.riv_y,
                                               self.n, self._nslope)

        # Floodplain sedimentation
        self.n, self.dn_fp = FP.dep_blanket(self._dy, self._dx, self._imax,
                                            self._jmax, self.current_SL,
                                            self.blanket_rate, self.n, self.riv_x,
                                            self.riv_y, self._ch_depth)

        # Wetland sedimentation
        ### no wetlands in first version of coupling to CEM ###
        self.n, self.dn_fp = FP.wetlands(self._dx, self._dy, self._imax, self._jmax,
                                         self.current_SL, self.WL_Z, self.WL_dist,
                                         self.n, self.riv_x, self.riv_y, self._x,
                                         self._y, self.dn_fp)

        # calculate sediment flux
        self.sed_flux = flux.calc_qs(self.nu, self.riv_x, self.riv_y, self.n,
                                     self._dx, self._dy, self.dt)

        # create a river profile array
        self.profile = prof.make_profile(self._dx, self._dy, self.n, self.riv_x,
                                         self.riv_y, self.profile)

        self.riv_mouth = [self.riv_x[-1], self.riv_y[-1]]

        self.k += 1
        self.time += self.dt

        # print "sediment flux = %f" % self.sed_flux

    def finalize(self):
        """Clean up & save avulsion file"""
        
        if self.savefiles == 1:
            np.savetxt('avulsions', self.avulsions, fmt='%i %i %i %.3f %.3f %.3f')
        pass

    def get_value_ref(self, var_name):
        return self._values[var_name]

    # def get_value(self, var_name):
    #     return self.get_value_ref(var_name).copy()

    def get_value_at_indices(self, var_name, indices):
        return self.get_value_ref(var_name).take(indices)

    def set_value(self, var_name, src):
        val = self.get_value_ref(var_name)
        val[:] = src

    def set_value_at_indices(self, var_name, src, indices):
        val = self.get_value_ref(var_name)
        val.flat[indices] = src

    def get_component_name(self):
        """Name of the component."""
        return self._name

    def get_input_var_names(self):
        """Get names of input variables."""
        return self._input_var_names

    def get_output_var_names(self):
        """Get names of output variables."""
        return self._output_var_names

    def get_current_time(self):
        return self._model.time/86400

    def get_time_step(self):
        """Time step of model."""
        return self._model.dt/86400

def main ():
    model = River_Module()
    model.initialize('input.yaml')

    while model.k < model.kmax:
        model.update()

        # save files
        if model.savefiles == 1:
            if model.k >= model.save_after:
                if model.k % model.savespacing == 0:
                    np.savetxt('elev_grid/elev_' + str(model.k*model.dt/86400 
                                - model.save_after) + '.out', model.n, fmt='%.6f')
                    np.savetxt('riv_course/riv_' + str(model.k*model.dt/86400
                                - model.save_after) + '.out',
                                zip(model.riv_x, model.riv_y), fmt='%i')
                    np.savetxt('profile/prof_' + str(model.k*model.dt/86400
                                - model.save_after) + '.out',
                                model.profile, fmt='%.6f')
                    np.savetxt('dn_fp/dn_fp_' + str(model.k*model.dt/86400
                                - model.save_after) + '.out',
                                model.dn_fp, fmt='%.6f')

    model.finalize()

if __name__ == '__main__':
    main()


