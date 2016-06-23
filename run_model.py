
#!/usr/bin/python
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import inspect, os
#from rafem.riverbmi import BmiRiverModule

N_DAYS = 5 * 365
Save_Daily_Timesteps = 1
Save_Yearly_Timesteps = 0
Save_Fluxes = 1
save_int = 1  # (in days)

def plot_coast(spacing, z):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    
    land = plt.cm.terrain(np.linspace(.4, 1, 128))
    ocean = plt.cm.ocean(np.linspace(.5, .8, 128))
    colors = np.vstack((ocean, land))
    m = LinearSegmentedColormap.from_list('land_ocean', colors)
    
    (x, y) = np.meshgrid(np.arange(z.shape[0]) * spacing[0],
                         np.arange(z.shape[1]) * spacing[1], indexing='ij')
    
    plt.pcolormesh(y * 1e-3, x * 1e-3, z, cmap=m, vmin=-50, vmax=50)
    
    plt.gca().set_aspect(1.) 
    plt.axis([0, 20, 0, 12])
    # plt.colorbar(orientation='horizontal').ax.set_xlabel('Elevation (m)')
    plt.xlabel('backwater lengths')
    plt.ylabel('backwater lengths')


from cmt.components import Cem, Rafem, Waves
cem = Cem()
raf = Rafem()
waves = Waves()

cem.setup('_run_cem', number_of_cols=200, number_of_rows=120, grid_spacing=100.)
raf.setup('_run_rafem', number_of_columns=200, number_of_rows=120, row_spacing=0.1,
          column_spacing=0.1, random_seed=623, rate_of_sea_level_rise=0.00,
          channel_discharge=10., upstream_elevation=5.)

cem.initialize('_run_cem/cem.txt')
raf.initialize('_run_rafem/input.yaml')
waves.initialize(None)

set(raf.get_output_var_names()) & set(cem.get_input_var_names())

z = raf.get_value('land_surface__elevation')
raf.set_value('land_surface__elevation', z)
cem.set_value('land_surface__elevation', z)

waves.set_value('sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_asymmetry_parameter', .5)
waves.set_value('sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_highness_parameter', .3)
cem.set_value("sea_surface_water_wave__height", 0.1)
cem.set_value("sea_surface_water_wave__period", 9.)
#cem.set_value("sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity", 0. * np.pi / 180.)

grid_id = cem.get_var_grid('land_surface__elevation')
spacing = cem.get_grid_spacing(grid_id)
shape = cem.get_grid_shape(grid_id)

z0 = raf.get_value('land_surface__elevation').reshape(shape)
riv_x = raf.get_value('channel_centerline__x_coordinate')/1000
riv_y = raf.get_value('channel_centerline__y_coordinate')/1000
# plot_coast(spacing, z0)
# plt.plot(riv_y,riv_x)
# plt.savefig('elev_initial.png')

qs = np.zeros_like(z0)
flux_array = np.zeros(2, dtype=np.float)

# np.random.seed(1991)

RIVER_WIDTH = dict(raf.parameters)['channel_width'] # Convert unit-width flux to flux
RHO_SED = 2650. # Used to convert volume flux to mass flux
TIME_STEP = raf.get_time_step()
Tcf = 1000000/86400

dx = (dict(raf.parameters)['row_spacing'])*1000.
slope = dict(raf.parameters)['delta_slope']
max_cell_h = dx*slope
channel_depth = dict(raf.parameters)['channel_depth']
max_rand = 0.0001

if Save_Daily_Timesteps or Save_Yearly_Timesteps:
# make directories to save run data
    if not os.path.exists("output_data_waves"):
        os.mkdir("output_data_waves")
    # if not os.path.exists("output_data_waves/elev_grid"):
    #     os.mkdir("output_data_waves/elev_grid")
    # if not os.path.exists("output_data_waves/riv_course"):
        # os.mkdir("output_data_waves/riv_course")
    # if not os.path.exists("output_data_waves/riv_profile"):
        # os.mkdir("output_data_waves/riv_profile")
    if not os.path.exists("output_data_waves/elev_figs"):
        os.mkdir("output_data_waves/elev_figs")
    if not os.path.exists("output_data_waves/prof_figs"):
        os.mkdir("output_data_waves/prof_figs")
    if not os.path.exists("output_data_waves/rel_elev"):
        os.mkdir("output_data_waves/rel_elev")

for time in np.arange(0, N_DAYS, TIME_STEP):

    raf.update_until(time)
    nyears = float(time/365.)

    sea_level = raf.get_value('sea_water_surface__elevation')
    
    # Get and set sediment flux at the river mouth
    raf_qs = raf.get_value('channel_exit_water_sediment~bedload__volume_flow_rate')
    
    y, x = raf.get_value('channel_exit__y_coordinate'), raf.get_value('channel_exit__x_coordinate')
    qs[int(y[0] / spacing[0]), int(x[0] / spacing[1])] = raf_qs[0] * RIVER_WIDTH * RHO_SED

    if Save_Fluxes:
        with open('fluxes.out','a') as file:
            file.write("%.2f %.5f \n" % (time, raf_qs[0] * RIVER_WIDTH * RHO_SED))

    # flux_array = np.vstack([flux_array, [(time, raf_qs[0] * RIVER_WIDTH * RHO_SED)]])
    
    cem.set_value('land_surface_water_sediment~bedload__mass_flow_rate', qs)

    # Get and set elevations from Rafem to CEM
    raf_z = (raf.get_value('land_surface__elevation') - sea_level).reshape(shape)
    riv_x = raf.get_value('channel_centerline__x_coordinate')/dx
    riv_y = raf.get_value('channel_centerline__y_coordinate')/dx
    riv_i = riv_x.astype(int)
    riv_j = riv_y.astype(int)
    prof_elev = raf_z[riv_j, riv_i]
    raf_z[riv_j, riv_i] += channel_depth

    # divide subaerial cells by max_cell_h to convert to percent full
    raf_z[raf_z > 0] /= max_cell_h

    # fix river elevations before passing
    mouth_cell_count = 0
    for k in reversed(xrange(riv_x.size)):
        if raf_z[riv_j[k], riv_i[k]] < 1:
            if mouth_cell_count < 1:
                mouth_cell_count += 1
            else:
                raf_z[riv_j[k], riv_i[k]] = 1

    # save grid passed to cem
    # np.savetxt('output_data_waves/cem_elev/cem_elev_'+str("%.3f" % nyears)+'.out',raf_z,fmt='%.5f')

    raf_z.reshape(shape[0]*shape[1],)
    cem.set_value('land_surface__elevation', raf_z)

    # update wave climate
    waves.update()
    angle = waves.get_value('sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity')

    cem.set_value('sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity', angle)
    cem.update_until(time)
    
    # Get and set elevations from CEM to Rafem
    cem_z = cem.get_value('land_surface__elevation').reshape(shape)
    cem_z[cem_z > 0] *= max_cell_h
    
    # reset river elevations back for Rafem
    # cem_z[[riv_j, riv_i] >= 0] -= channel_depth
    if cem_z[riv_j[-1],riv_i[-1]] > 0:
        cem_z[riv_j[:-1],riv_i[:-1]] = prof_elev[:-1]
        cem_z[riv_j[-1],riv_i[-1]] -= channel_depth
    
    else:
        cem_z[riv_j[:-2],riv_i[:-2]] = prof_elev[:-2]
        cem_z[riv_j[-2],riv_i[-2]] -= channel_depth

    cem_z.reshape(shape[0]*shape[1],)
    raf.set_value('land_surface__elevation', cem_z + sea_level)
    
    qs.fill(0)

    if time % save_int == 0:

        print('time = %.3f days' % time)
        # save outputs
        z = raf.get_value('land_surface__elevation').reshape(shape)
        rel_z = z - sea_level
        x = raf.get_value('channel_centerline__x_coordinate')
        y = raf.get_value('channel_centerline__y_coordinate')
        prof = raf.get_value('channel_centerline__elevation')
        real_prof = rel_z[(y/dx).astype(int), (x/dx).astype(int)]
        river_x = x/1000
        river_y = y/1000
        riv_left = z[y.astype(int)/100, (x.astype(int)/100) + 1]
        riv_right = z[y.astype(int)/100, (x.astype(int)/100) - 1]
        riv_left[riv_left < sea_level] = sea_level
        riv_right[riv_right < sea_level] = sea_level
        Tcf_time = time/Tcf



        ### SAVE DAILY TIMESTEPS ###
        ##########################################################################################
        if Save_Daily_Timesteps == 1:

            # np.savetxt('output_data_waves/elev_grid/elev_'+str("%.3f" % nyears)+'.out',z,fmt='%.5f')
            np.savetxt('output_data_waves/rel_elev/rel_elev_'+str("%i" % time)+'.out',rel_z,fmt='%.5f')
            # np.savetxt('output_data_waves/riv_course/riv_'+str("%i" % time)+'.out',zip(x,y),fmt='%i')
            # np.savetxt('output_data_waves/riv_profile/prof_'+str("%i" % time)+'.out',real_prof,fmt='%.5f')

            # save figures
            f = plt.figure()
            plot_coast(spacing, z - sea_level)
            plt.plot(river_x, river_y, LineWidth=2.5)
            plt.title('time = '+str("%.3f" % Tcf_time)+' Tcf')
            plt.savefig('output_data_waves/elev_figs/elev_fig_'+str(int(time))+'.png')
            plt.close(f)

            p = plt.figure()
            PL, = plt.plot(riv_left, color=[0, 0.2, 0], linewidth=2)
            SL = plt.hlines(sea_level, 0, 120, color='c', linewidth=1)
            PP, = plt.plot(prof, color='b', linewidth=2)
            PR, = plt.plot(riv_right, color=[0, 0.8, 0], linewidth=2, linestyle='dashed')
            plt.legend((PL, PR, PP, SL),('adjacent floodplain left side',
                        'adjacent floodplain right side', 'river bed', 
                        'sea level'), fontsize=14, loc=('upper right'))
            plt.axis([0, 120, -10, 20])
            plt.title('time = '+str("%.3f" % Tcf_time)+' Tcf')
            plt.xlabel('river cells')
            plt.ylabel('channel depths')
            plt.savefig('output_data_waves/prof_figs/prof_fig_'+str(int(time))+'.png')
            plt.close(p)
        ##########################################################################################
        

        ### SAVE YEARLY TIMESTEPS ###
        ##########################################################################################
        if Save_Yearly_Timesteps == 1:
            # np.savetxt('output_data_waves/elev_grid/elev_'+str(time/save_int)+'.out',z,fmt='%.5f')
            np.savetxt('output_data_waves/rel_elev/rel_elev_'+str(time/save_int)+'.out',rel_z,fmt='%.5f')
            np.savetxt('output_data_waves/riv_course/riv_'+str(time/save_int)+'.out',zip(x,y),fmt='%i')
            np.savetxt('output_data_waves/riv_profile/prof_'+str(time/save_int)+'.out',real_prof,fmt='%.5f')

            # save figures
            f = plt.figure()
            plot_coast(spacing, z - sea_level)
            plt.plot(river_y,river_x, linewidth=2)
            plt.title('time = '+str(time/save_int)+' years, sea level = '+str("%.3f" % sea_level)+' m')
            plt.savefig('output_data_waves/elev_figs/elev_fig_'+str(time/save_int)+'.png')
            plt.close(f)

            p = plt.figure()
            plt.plot(prof)
            plt.axis([0, 200, -20, 120])
            plt.title('time = '+str(time/save_int)+' years')
            plt.savefig('output_data_waves/prof_figs/prof_fig_'+str(time/save_int)+'.png')
            plt.close(p)
        ##########################################################################################

# np.savetxt('output_data_waves/fluxes.out',flux_array,fmt='%i,%.5f')
