
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import inspect, os
import pdb
#from rafem.riverbmi import BmiRiverModule

def plot_coast(spacing, z):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    
    land = plt.cm.terrain(np.linspace(.4, 1, 128))
    ocean = plt.cm.ocean(np.linspace(.5, .8, 128))
    colors = np.vstack((ocean, land))
    m = LinearSegmentedColormap.from_list('land_ocean', colors)
    
    (y, x) = np.meshgrid(np.arange(z.shape[0]) * spacing[0],
                         np.arange(z.shape[1]) * spacing[1], indexing='ij')
    
    plt.pcolormesh(y * 1e-3, x * 1e-3, z, cmap=m, vmin=-100, vmax=100)
    
    plt.gca().set_aspect(1.) 
    plt.colorbar(orientation='horizontal').ax.set_xlabel('Elevation (m)')
    plt.xlabel('Cross-shore (km)')
    plt.ylabel('Alongshore (km)')


from cmt.components import Cem, Rafem, Waves
cem = Cem()
raf = Rafem()
waves = Waves()

cem.setup('_run_cem', number_of_cols=100, number_of_rows=300, grid_spacing=5000.)
raf.setup('_run_rafem', number_of_columns=100, number_of_rows=300, row_spacing=5.,
          column_spacing=5., rate_of_sea_level_rise=0.00, channel_discharge=10000.)

cem.initialize('_run_cem/cem.txt')
raf.initialize('_run_rafem/input.yaml')
waves.initialize(None)

set(raf.get_output_var_names()) & set(cem.get_input_var_names())

z = raf.get_value('land_surface__elevation')
z = z - 0.05
raf.set_value('land_surface__elevation', z)
cem.set_value('land_surface__elevation', z)

waves.set_value('sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_asymmetry_parameter', .5)
waves.set_value('sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_highness_parameter', .2)
cem.set_value("sea_surface_water_wave__height", 0.05)
cem.set_value("sea_surface_water_wave__period", 9.)
#cem.set_value("sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity", 0. * np.pi / 180.)

grid_id = cem.get_var_grid('land_surface__elevation')
spacing = cem.get_grid_spacing(grid_id)
shape = cem.get_grid_shape(grid_id)

z0 = raf.get_value('land_surface__elevation').reshape(shape)
riv_x = raf.get_value('channel_centerline__x_coordinate')/1000
riv_y = raf.get_value('channel_centerline__y_coordinate')/1000
plot_coast(spacing, z0)
plt.plot(riv_y,riv_x)
# plt.savefig('elev_initial.png')

qs = np.zeros_like(z0)
flux_array = np.zeros(2, dtype=np.float)

RIVER_WIDTH = dict(raf.parameters)['channel_width'] # Convert unit-width flux to flux
RHO_SED = 2650. # Used to convert volume flux to mass flux
N_DAYS = 100 * 365
TIME_STEP = int(raf.get_time_step())
save_int = 365

dx = (dict(raf.parameters)['row_spacing'])*1000.
slope = dict(raf.parameters)['delta_slope']
max_cell_h = dx*slope
channel_depth = dict(raf.parameters)['channel_depth']
max_rand = 0.000001

riv_length = len(riv_x)

# make directories to save run data
import os
if not os.path.exists("output_data_waves"):
    os.mkdir("output_data_waves")
if not os.path.exists("output_data_waves/elev_grid"):
    os.mkdir("output_data_waves/elev_grid")
if not os.path.exists("output_data_waves/riv_course"):
    os.mkdir("output_data_waves/riv_course")
if not os.path.exists("output_data_waves/riv_profile"):
    os.mkdir("output_data_waves/riv_profile")
if not os.path.exists("output_data_waves/elev_figs"):
    os.mkdir("output_data_waves/elev_figs")
if not os.path.exists("output_data_waves/prof_figs"):
    os.mkdir("output_data_waves/prof_figs")
if not os.path.exists("output_data_waves/rel_elev"):
    os.mkdir("output_data_waves/rel_elev")
#if not os.path.exists("output_data_waves/shoreline_pos"):
#    os.mkdir("output_data_waves/shoreline_pos")
#if not os.path.exists("output_data_waves/lowest_cell"):
#    os.mkdir("output_data_waves/lowest_cell")

def lowest_adj_cell(n, sub):
    i,j = sub

    if j == 0 and i == 0:
        di, dj = np.array([1, 1, 0]), np.array([0, 1, 1])
    elif j == 0 and i == n.shape[0] - 1:
        di, dj = np.array([-1, -1, 0]), np.array([0, 1, 1])
    elif j == n.shape[1] - 1 and i == 0:
        di, dj = np.array([0, 1, 1]), np.array([-1, -1, 0])
    elif j == n.shape[1] - 1 and i == n.shape[0] - 1:
        di, dj = np.array([0, -1, -1]), np.array([-1, -1, 0])
    elif j == n.shape[1] - 1:
        di, dj  = np.array([-1, -1, 0, 1, 1]), np.array([0, -1, -1, -1, 0])
    elif j == 0:
        di, dj  = np.array([-1, -1, 0, 1, 1]), np.array([0, 1, 1, 1, 0])
    elif i == n.shape[0] - 1:
        di, dj = np.array([0, -1, -1, -1, 0]), np.array([-1, -1, 0, 1, 1])
    elif i == 0:
        di, dj = np.array([0, 1, 1, 1, 0]), np.array([-1, -1, 0, 1, 1])
    else:
        di, dj = np.array([0, -1, -1, -1, 0, 1, 1, 1]),  np.array([-1, -1, 0, 1, 1, 1, 0, -1])

    lowest = np.amin(n[i + di, j + dj])

    return lowest

for time in xrange(0, N_DAYS, TIME_STEP):
    # if time % 5000 == 0:
    #     print 'running rafem until {time}'.format(time=time)

    raf.update_until(time)

    sea_level = raf.get_value('sea_water_surface__elevation')
    
    # Get and set sediment flux at the river mouth
    raf_qs = raf.get_value('channel_exit_water_sediment~bedload__volume_flow_rate')
    
    y, x = raf.get_value('channel_exit__y_coordinate'), raf.get_value('channel_exit__x_coordinate')
    qs[int(y[0] / spacing[0]), int(x[0] / spacing[1])] = raf_qs[0] * RIVER_WIDTH * RHO_SED
    
    flux_array = np.vstack([flux_array, [(time, raf_qs[0] * RIVER_WIDTH * RHO_SED)]])
    
    cem.set_value('land_surface_water_sediment~bedload__mass_flow_rate', qs)
    
    # np.savetxt('qs.out',qs,'%.3f')

    # Get and set elevations from Rafem to CEM
    
    raf_z = (raf.get_value('land_surface__elevation') - sea_level).reshape(shape)
    # np.savetxt('raf_z_orig.out',raf_z,'%.5f')
    riv_x = raf.get_value('channel_centerline__x_coordinate')/dx
    riv_y = raf.get_value('channel_centerline__y_coordinate')/dx
    riv_i = riv_x.astype(int)
    riv_j = riv_y.astype(int)
    # if len(riv_x) != riv_length:
    #     pdb.set_trace()
    prof_elev = raf_z[riv_j, riv_i]
    raf_z[riv_j, riv_i] += channel_depth
    # raf_z[riv_j[:-1], riv_i[:-1]] += 2*channel_depth
    # raf_z[riv_j[-1], riv_i[-1]] += channel_depth

    #np.savetxt('extracted_prof.out',raf.get_value('channel_centerline__elevation'),'%.3f')
    #np.savetxt('prof_before_cem.out',prof_elev,'%.3f')
    #np.savetxt('raf_z_ch_added.out',raf_z,'%.5f')
    
    # # find shoreline cells and fix elevations of subaerial cells
    # shore_flag = np.zeros(raf_z.shape, dtype=np.int)
    # lowest_cell = np.zeros(raf_z.shape)
    # for i in xrange(raf_z.shape[0]):
    #     for j in xrange(raf_z.shape[1]):
    #         #lowest_elev = lowest_adj_cell(raf_z, (i,j))
    #         #lowest_cell[i,j] = lowest_elev
    #         if 0 < raf_z[i,j] <= max_cell_h:
    #             if lowest_adj_cell(raf_z, (i,j)) <= 0:
    #                 shore_flag[i,j] = 1
    #             else:
    #                 raf_z[i,j] = max_cell_h + np.random.rand()*max_rand
    #                 shore_flag[i,j] = 2

    # divide subaerial cells by max_cell_h to convert to percent full
    raf_z[raf_z > 0] /= max_cell_h
    
    #np.savetxt('raf_z_divided.out',raf_z,'%.5f')
    
    #np.savetxt('output_data_waves/shoreline_pos/shoreline'+str(time)+'.out',shore_flag,'%i')
    #np.savetxt('output_data_waves/lowest_cell/low_cell'+str(time)+'.out',lowest_cell,'%.3f')

    # fix river elevations before passing
    mouth_cell_count = 0
    for k in reversed(xrange(riv_x.size)):
        if raf_z[riv_j[k], riv_i[k]] < 1:
            if mouth_cell_count < 1:
                mouth_cell_count += 1
            else:
                raf_z[riv_j[k], riv_i[k]] = 1
    
    #np.savetxt('river_depitted.out',raf_z[riv_j,riv_i],'%.3f')
    #np.savetxt('raf_z_ready_to_pass.out',raf_z,'%.3f')
    
    # if len(riv_x) != riv_length:
    #     pdb.set_trace()
    #     riv_length = len(riv_x)

    raf_z.reshape(shape[0]*shape[1],)
    cem.set_value('land_surface__elevation', raf_z)

    # update wave climate
    waves.update()
    angle = waves.get_value('sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity')

    cem.set_value('sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity', angle)
    cem.update_until(time)
    
    # Get and set elevations from CEM to Rafem
    cem_z = cem.get_value('land_surface__elevation').reshape(shape)
    #np.savetxt('cem_z_raw.out',cem_z,'%.3f')
    
    cem_z[cem_z > 0] *= max_cell_h
    
    #np.savetxt('cem_z_divided.out',cem_z,'%.3f')
    
    # cem_z[[riv_j, riv_i] >= 0] -= channel_depth
    
    if cem_z[riv_j[-1],riv_i[-1]] > 0:
        cem_z[riv_j[:-1],riv_i[:-1]] = prof_elev[:-1]
        cem_z[riv_j[-1],riv_i[-1]] -= channel_depth
    
    else:
        cem_z[riv_j[:-2],riv_i[:-2]] = prof_elev[:-2]
        cem_z[riv_j[-2],riv_i[-2]] -= channel_depth

    #np.savetxt('reset_cem_prof.out',cem_z[riv_j,riv_i],'%.3f')
    #np.savetxt('cem_to_rafem.out',cem_z,'%.3f')
    
    cem_z.reshape(shape[0]*shape[1],)
    raf.set_value('land_surface__elevation', cem_z + sea_level)
    
    qs.fill(0)

    if time % save_int == 0:
        # save outputs
        z = raf.get_value('land_surface__elevation').reshape(shape)
        rel_z = z - sea_level
        x = raf.get_value('channel_centerline__x_coordinate')
        y = raf.get_value('channel_centerline__y_coordinate')
        prof = raf.get_value('channel_centerline__elevation')
        real_prof = rel_z[(y/dx).astype(int), (x/dx).astype(int)]

        np.savetxt('output_data_waves/elev_grid/elev_'+str(time/save_int)+'.out',z,fmt='%.5f')
        np.savetxt('output_data_waves/rel_elev/rel_elev_'+str(time/save_int)+'.out',rel_z,fmt='%.5f')
        np.savetxt('output_data_waves/riv_course/riv_'+str(time/save_int)+'.out',zip(x,y),fmt='%i')
        np.savetxt('output_data_waves/riv_profile/prof_'+str(time/save_int)+'.out',real_prof,fmt='%.5f')

        river_x = x/1000
        river_y = y/1000

        # save figures
        f = plt.figure()
        plot_coast(spacing, z - sea_level)
        plt.plot(river_y,river_x)
        plt.title('time = '+str(time/save_int)+' years')
        plt.savefig('output_data_waves/elev_figs/elev_fig_'+str(time/save_int)+'.png')
        plt.close(f)

        p = plt.figure()
        plt.plot(prof)
        plt.axis([0, 300, -20, 120])
        plt.title('time = '+str(time/save_int)+' years')
        plt.savefig('output_data_waves/prof_figs/prof_fig_'+str(time/save_int)+'.png')
        plt.close(p)
    
#    if cem_z[riv_j[-2],riv_i[-2]] < cem_z[riv_j[-1],riv_i[-1]]:
#        break
        
np.savetxt('output_data_waves/fluxes.out',flux_array,fmt='%i,%.5f')
