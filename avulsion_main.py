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

def initialize(fname):
    """ initialize the avulsion module """

    params = read_params_from_file('input.yaml')

    # Spatial parameters
    length, width = params['shape']
    dx_km, dy_km = params['spacing']
    L = length * 1000       # convert to meters
    W = width * 1000
    dx = dx_km * 1000
    dy = dy_km * 1000

    imax = L/dx + 1
    jmax = W/dy + 1
    x = np.zeros((imax, jmax))   # longitudinal space
    y = np.zeros((imax, jmax))   # transverse space
    n = np.zeros((imax, jmax))   # eta, elevation
    dn_rc = np.zeros((imax))       # change in elevation along river course
    dn_fp = np.zeros((imax, jmax))     # change in elevation due to floodplain dep
    riv_x = [0]             # defines first x river locations
    riv_y = [W/2]          # defines first y river locations
    profile = np.zeros((imax))  # elevation profile of river course
    avulsions = [(0, 0, 0, 0, 0, 0)]    # initializes timestep/avulsions array

    # Time parameters
    dt = (params['dt_day'] *60*60*24)     # convert timestep to seconds
    time_max_s = (params['time_max'] * 31536000)  # length of model run in seconds
    spinup_s = (params['spinup'] * 31536000)  # length of spinup in seconds
    kmax = spinup_s/dt + time_max_s/dt + 1  # max number of timesteps
    save_after = spinup_s/dt        # save files after this point
    time = 0.
    k = 0

    # Sea level and subsidence parameters
    SL = [params['Initial_SL']]                   # initializes SL array
    SLRR = (params['SLRR_m'] / 31536000) * dt  # sea level rise rate in m/s per timestep
    IRR = (params['IRR_m'] / 31536000) * dt    # inlet rise rate in m/s per timestep

    # River parameters
    init_cut = params['init_cut_frac'] * params['ch_depth']

    # Floodplain and wetland characteristics
    blanket_rate = (params['blanket_rate_m'] / 31536000) * dt    # blanket deposition in m/s
    splay_dep = (params['splay_dep_m'] / 31536000) * dt       # splay deposition in m/s

    # Initialize elevation grid
    for i in range(imax):
        for j in range(jmax):
            x[i][j] = i * dx
            y[i][j] = j * dy
            n[i][j] = params['n0'] - (params['nslope'] * float(x[i][j]) \
                      + params['max_rand'] * np.random.rand())
            j += 1
        i += 1

    params['L'] = L
    params['W'] = W
    params['dx'] = dx
    params['dy'] = dy
    params['imax'] = imax
    params['jmax'] = jmax
    params['x'] = x
    params['y'] = y
    params['n'] = n
    params['dn_rc'] = dn_rc
    params['dn_fp'] = dn_fp
    params['riv_x'] = riv_x
    params['riv_y'] = riv_y
    params['profile'] = profile
    params['SL'] = SL
    params['avulsions'] = avulsions
    params['dt'] = dt
    params['time_max_s'] = time_max_s
    params['spinup_s'] = spinup_s
    params['kmax'] = kmax
    params['save_after'] = save_after
    params['time'] = time
    params['k'] = k
    params['SLRR'] = SLRR
    params['IRR'] = IRR
    params['init_cut'] = init_cut
    params['blanket_rate'] = blanket_rate
    params['splay_dep'] = splay_dep

    # Determine initial river course
    riv_x, riv_y = steep_desc.find_course(params['dx'], params['dy'], params['imax'],
                   params['jmax'], params['n'], params['riv_x'], params['riv_y'])
    params['riv_x'] = riv_x
    params['riv_y'] = riv_y

    # downcut into new river course by amount determined by init_cut
    n = downcut.cut_init(params['dx'], params['dy'], params['riv_x'], params['riv_y'],
        params['n'], params['init_cut'])
    params['n'] = n

    # smooth initial river course elevations using linear diffusion equation
    n, dn_rc = diffuse.smooth_rc(params['dx'], params['dy'], params['nu'], params['dt'],
               params['riv_x'], params['riv_y'], params['n'], params['nslope'])
    params['n'] = n
    params['dn_rc'] = dn_rc

    # Determine initial river profile
    profile = prof.make_profile(params['dx'], params['dy'], params['n'], params['riv_x'],
              params['riv_y'], params['profile'])
    params['profile'] = profile

    # make directories and save initial condition files
    if params['savefiles'] == 1:
        # os.mkdir("run" + str(run_num) + "_out")
        os.mkdir("elev_grid")
        os.mkdir("riv_course")
        os.mkdir("profile")
        os.mkdir("dn_fp")
    #   saves initial conditions
    #    np.savetxt('elev_grid/elev_0.out', n, fmt='%f')
    #    np.savetxt('riv_course/riv_0.out', zip(riv_x, riv_y), fmt='%i')
    #    np.savetxt('profile/prof_0.out', profile, fmt='%f')

    return params

def update(params):
    """ Update avulsion model one time step. """

    # begin time loop and main program
    # for k in range(kmax):

    # determine current sea level
    SL = params['SL'] + [params['k'] * params['SLRR']]
    current_SL = SL[-1]
    params['SL'] = SL
    params['current_SL'] = current_SL

    ### future work: SLRR can be a vector to change rates ###

    # determine if there is an avulsion & find new path if so
    riv_x, riv_y, loc, SEL, SER, n, dn_fp, avulsion_type, length_new_sum, \
        length_old = avulse.find_avulsion(params['dx'], params['dy'], 
                        params['imax'], params['jmax'], params['riv_x'],
                        params['riv_y'], params['n'], params['super_ratio'],
                        params['current_SL'], params['ch_depth'],
                        params['short_path'], params['dn_fp'], params['splay_type'],
                        params['splay_dep'])
    params['riv_x'] = riv_x
    params['riv_y'] = riv_y
    params['SEL'] = SEL
    params['SER'] = SER
    params['n'] = n
    params['dn_fp'] = dn_fp

    # save timestep and avulsion location if there was one
    if len(loc) != 0:
        avulsions = (params['avulsions'] + [(params['k']*params['dt']/86400,
                    loc[-1], avulsion_type, length_old, length_new_sum,
                    current_SL)])
        params['avulsions'] = avulsions
    
    # raise first two rows by inlet rise rate (subsidence)
    params['n'][0][:] = params['n'][0][:] + (params['IRR'])
    params['n'][1][:] = params['n'][1][:] + (params['IRR'])

    # change elevations according to sea level rise (SLRR)
    n, rc_flag = SLR.elev_change(params['imax'], params['jmax'],
                    params['current_SL'], params['n'], params['riv_x'],
                    params['riv_y'], params['ch_depth'], params['dx'],
                    params['dy'])
    params['n'] = n
    params['rc_flag'] = rc_flag

    # smooth river course elevations using linear diffusion equation
    n, dn_rc = diffuse.smooth_rc(params['dx'], params['dy'], params['nu'],
                params['dt'], params['riv_x'], params['riv_y'], params['n'],
                params['nslope'])
    params['n'] = n
    params['dn_rc'] = dn_rc

    # Floodplain sedimentation
    n, dn_fp = FP.dep_blanket(params['dy'], params['dx'], params['imax'],
                params['jmax'], params['current_SL'], params['blanket_rate'],
                params['n'], params['riv_x'], params['riv_y'], params['ch_depth'])
    params['n'] = n
    params['dn_fp'] = dn_fp
    
    # Wetland sedimentation
    ### no wetlands in first version of coupling to CEM ###
    n, dn_fp = FP.wetlands(params['dx'], params['dy'], params['imax'], 
                params['jmax'], params['current_SL'], params['WL_Z'], 
                params['WL_dist'], params['n'], params['riv_x'], params['riv_y'],
                params['x'], params['y'], params['dn_fp'])
    params['n'] = n
    params['dn_fp'] = dn_fp

    # calculate sediment flux
    sed_flux = flux.calc_qs(params['nu'], params['riv_x'], params['riv_y'], params['n'],
                params['dx'], params['dy'], params['dt'])
    params['sed_flux'] = sed_flux

    # create a river profile array
    profile = prof.make_profile(params['dx'], params['dy'], params['n'], params['riv_x'],
              params['riv_y'], params['profile'])
    params['profile'] = profile

    # save files
    if params['savefiles'] == 1:
        if params['k'] >= params['save_after']:
            if params['k'] % params['savespacing'] == 0:
                np.savetxt('elev_grid/elev_' + str(params['k']*params['dt']/86400 
                            - (params['save_after'])) + '.out', params['n'], fmt='%.6f')
                np.savetxt('riv_course/riv_' + str(params['k']*params['dt']/86400
                            - (params['save_after'])) + '.out',
                            zip(params['riv_x'], params['riv_y']), fmt='%i')
                np.savetxt('profile/prof_' + str(params['k']*params['dt']/86400
                            - (params['save_after'])) + '.out',
                            params['profile'], fmt='%.6f')
                np.savetxt('dn_fp/dn_fp_' + str(params['k']*params['dt']/86400
                            - (params['save_after'])) + '.out',
                            params['dn_fp'], fmt='%.6f')
    if params['savefiles'] == 1:
        np.savetxt('avulsions', params['avulsions'], fmt='%.3f')

    params['k'] += 1
    params['time'] += params['dt']

def finalize():
    """Finalize the avulsion model."""
    pass

def main ():

    params = initialize('input.yaml')

    while params['k'] < params['kmax']:
        update(params)

    # print "sediment flux = %f" % params['sed_flux']

    finalize()

if __name__ == '__main__':
    main()


