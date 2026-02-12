import numpy as np
# import h5py as h5
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
# from scipy.ndimage.filters import gaussian_filter
# import matplotlib.patches as pat
# plt.style.use('dark_background')
# plt.style.use('default')
import warnings; warnings.simplefilter('ignore')
cmap=cm.viridis
from scipy.optimize import curve_fit
# from matplotlib.lines import Line2D
# import matplotlib.patches as mpatches
# from matplotlib import gridspec
# from scipy import stats
# from astropy.cosmology import FlatLambdaCDM
# import astropy.units as u
# import math
# from copy import copy
# from matplotlib.patches import Circle
# from matplotlib.patches import Rectangle
# from collections import Counter
# # from sklearn.decomposition import nearest_lparticle
# from sklearn.preprocessing import StandardScaler
# import numpy.linalg as lina
# import os.path
# from scipy.fft import fft, ifftn, fftn, fftfreq
# import struct
# import ctypes
# import random
import os
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler
# import nbodykit.lab as nb
# from nbodykit.source.catalog import ArrayCatalog
# from nbodykit.algorithms.fftcorr import FFTCorr as fftcor
# sys.path.append('/cosma7/data/dp004/dc-wang4/Nexus/NEXUS_1.1/python/python_import/')
# import MMF
# from sklearn.neighbors import KernelDensity
# from matplotlib import rcParams
from matplotlib import rc
# rc('text', usetex=True)
# rc('font',**{'family':'serif'})
rc('font',**{'family':'serif','serif':['Times']})
# from matplotlib.ticker import FormatStrFormatter
import pandas as pd
# from scipy.stats.stats import spearmanr  
# from scipy import signal
# from scipy import interpolate as interp
# from scipy.signal import savgol_filter
# from scipy.stats.stats import pearsonr  
import time
sys.path.append('/cosma/home/durham/dc-wang4/Python_functions/')
from General_functions import *
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')


#t1=time.clock()

L = int(sys.argv[1])
HaloID = int(sys.argv[2])
kfs_list = [7, 10, 12, 20, 24, 30, 35, 50]
ik = int(sys.argv[3])
kfs = kfs_list[ik]
halor = 4
T_min = int(sys.argv[4])
T_max = int(sys.argv[5])
Nh = int(sys.argv[6])
Nf = int(sys.argv[7])
Nfile = int(sys.argv[8])
# side_length = int(sys.argv[8])

print('Level: ', L)
print('HaloID: ', HaloID)
print('Nh: ', Nh)
print('Nf: ', Nf)
print('T min: ', T_min)
print('T max: ', T_max)
print('kfs: ', kfs)

if HaloID == 108:
    side_length = 1.317
elif HaloID == 150:
    side_length = 1.456
elif HaloID == 225:
    side_length = 0.956
elif HaloID == 206:
    side_length = 1.100
elif HaloID == 260:
    side_length = 1.057
elif HaloID == 305:
    side_length = 1.015
elif HaloID == 337:
    side_length = 0.783
elif HaloID == 443:
    side_length = 1.100
    
print('side_length: ', side_length)

# Indexer_Halo = {}
# Halo_Dist = {}
# Dist_Axis = {}
# Halo_Density = {}
# Par_SubIndex = {}

# for T_min_temp in range(0, T_max+1):
#     for T_max_temp in range(200, T_min_temp, -1):
#         label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min_temp, T_max_temp)
#         path = '/cosma7/data/dp004/dc-wang4/VVV_data/Halo_Data/AllSubRef_{}.npy'.format(label)
#         try:
#             AllHalo_data = np.load(path, allow_pickle = True)
            
#             Indexer_Halo.update(AllHalo_data[0])
#             Halo_Dist.update(AllHalo_data[1])
#             Dist_Axis.update(AllHalo_data[2])
#             Halo_Density.update(AllHalo_data[3])
#             Par_SubIndex.update(AllHalo_data[4])
#             print(path)
#         except:
#             pass

Hconst = 0.677

Par_Coor = {}
Par_Mass = {}
Par_ID = {}
Par_Vel = {}
Halo_Mass_200 = {}
SubHalo_Mass = {}
Halo_FirSub = {}
Halo_NSub = {}
SubHalo_GrN = {}
Halo_Offset = {}
SubHalo_Offset = {}
Halo_ParN = {}
SubHalo_ParN = {}

Halo_Mass = {}
Halo_Coor = {}
Halo_R_200 = {}
Subhalo_Coor = {}
Halo_Vel = {}

Subhalo_HalfmassRad = {}
Subhalo_IDMostbound = {}
Subhalo_Spin = {}
Subhalo_Vel = {}
Subhalo_VelDisp = {}
Subhalo_Vmax = {}
Subhalo_VmaxRad = {}

Redshift = {}
Box_size = {}
Par_LinearDensity = {}
# Subhalo_FirstProg = {}

print('finish initializing')

linearDensity = True
for T in [-1]:
    if Nh > 15:
        path = '/cosma7/data/dp004/dc-wang4/IC_gen/apply_ps/output/' + \
        'L{}_HaloID{}Nh{}Nf{}Nc8r{}k{}/output_2lpt/ini_pos'.format(L, HaloID, Nh, Nf, halor, kfs)
    elif Nh == 15:
        path = '/cosma7/data/dp004/dc-wang4/IC_gen/apply_ps/output/' + \
        'L{}_HaloID{}Nh{}Nf{}Nc8r{}k{}/output/ini_pos'.format(L, HaloID, Nh, Nf, halor, kfs)
    
    print(path)
    File = Read_file_ini(path, 256, linearDensity)
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
    
    Par_Coor[label] = File[0]/Hconst
    Par_Mass[label] = File[1]*10**10/Hconst
    Par_ID[label] = File[2]
    # Par_Vel[label] = File[3]/Hconst

    if linearDensity is True:
        Par_LinearDensity[label] = File[4]
        Redshift[label] = File[5]
    else:
        Redshift[label] = File[4]

side_length = side_length/Hconst

rfs = np.pi/kfs/Hconst
max_npeak = int((side_length/2/rfs)**3)
print('maximum number of peak:', max_npeak)


NbinH = 70
Nrand = int(1e3)
dens_factorH = 100
Nbin = 100
Ndistbin = 100
Nrhobin = 100
N_slice = 0

label0 = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, -1)
Par_x0 = Par_Coor[label0][:, 0]
Par_y0 = Par_Coor[label0][:, 1]
Par_z0 = Par_Coor[label0][:, 2]
Par_m0 = Par_Mass[label0]
Par_id0 = Par_ID[label0]
Par_linrho = Par_LinearDensity[label0]
index_parent = pd.Index(Par_id0)
        
for T in range(T_min, T_max+1):
    
    if Nfile == 6:
        if HaloID == 150 and Nh == 15:
            name = 'L{}_HaloID{}Nh{}Nf{}Nc8r{}k{}/'.format(L, HaloID, Nh, Nf, halor, kfs)
        else:
            name = 'L{}_HaloID{}Nh{}Nf{}Nc8r{}k{}/output/'.format(L, HaloID, Nh, Nf, halor, kfs)
        Path = '/cosma6/data/dp004/dc-wang4/' + name
    elif Nfile == 7:
        name = 'L{}_HaloID{}Nh{}Nf{}Nc8r{}k{}/output/'.format(L, HaloID, Nh, Nf, halor, kfs)
        Path = '/cosma7/data/dp004/dc-wang4/gadget_runs/' + name
    
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
    print(Path, label)
    try:
        File = Read_file(L, T, 16, Path)

        Par_Coor[label] = File[0]/Hconst
        Par_Mass[label] = File[1]*10**10/Hconst
        Par_ID[label] = File[2]
        # Par_Vel[label] = File[3]/Hconst
        Halo_Mass_200[label] = File[4]*10**10/Hconst
        # SubHalo_Mass[label] = File[5]*10**10/Hconst
        # Halo_FirSub[label] = File[6]
        Halo_NSub[label] = File[7]
        # SubHalo_GrN[label] = File[8]
        Halo_Offset[label] = File[9]
        Halo_ParN[label] = File[10]
        # SubHalo_Offset[label] = File[11]
        # SubHalo_ParN[label] = File[12]
        Halo_Mass[label] = File[13]*10**10/Hconst
        Halo_Coor[label] = File[14]/Hconst
        # Subhalo_Coor[label] = File[15]/Hconst
        Halo_R_200[label] = File[16]/Hconst
        Redshift[label] = File[17]
        Box_size[label] = File[18]/Hconst
        # Halo_Vel[label] = File[19]/Hconst
        # Subhalo_Spin[label] = File[20]
        
    except:
        print('T_max = ', T-1)
        T_max = T-1
        break
# 
print('finished reading', T_min, T_max)
if T_max < T_min:
    print('no file found')
    exit()

SubProg_Index = {}
SubProg_Mass = {}
Halo_Index = {}
label = 'L{}HaloID{}k{}Nh{}'.format(L, HaloID, kfs, Nh)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/Prog_Index_{}.npy'.format(label)
try:
    MergeTree_data = np.load(path, allow_pickle = True)
    SubProg_Index.update(MergeTree_data[0])
    SubProg_Mass.update(MergeTree_data[1])
    Halo_Index.update(MergeTree_data[2])
    print(path)
except:
    pass

NPeak_Test = {}
Peak_List = {}
Grid_Sizes = {}
label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_min)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Peak_Data/Identifier_{}.npy'.format(label)
try:
    AllPeak_data = np.load(path, allow_pickle = True)
    NPeak_Test.update(AllPeak_data[0])
    Peak_List.update(AllPeak_data[1])
    Grid_Sizes.update(AllPeak_data[2])
    print(path)
except:
    pass

Peaks_Properties = {}
CoorAllHaloPar0 = {}
LinRhoAllHaloPar = {}
Peak_IDs = {}
Semi_Selection = {}

for T in range(T_min, T_max+1):  
    plt.clf()
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
    print(label)

    if len(Halo_Mass_200[label]) == 0:
        print('No Halo at ', label)
        continue

    # T_index_selection = sub_prog_T == T
    # if np.sum(T_index_selection) == 0:
    #     print('No progenitor at ', label)
    #     continue
    try:
        halo_index = int(Halo_Index[label])
    except:
        print('No progenitor at ', label)
        continue

    Par_x = Par_Coor[label][:, 0]
    Par_y = Par_Coor[label][:, 1]
    Par_z = Par_Coor[label][:, 2]
    Par_id = Par_ID[label]

    Par_id_unique_list = np.unique(Par_id)
    # print(len(Par_id_unique_list), len(Par_id), 
    #       len(Par_id_unique_list)==len(Par_id))
    # exit()

    nhalo = len(Halo_Mass_200[label])
    if nhalo == 0:
        T_min = T+1
        print('no halohalo_nsub in: ', label)
        continue
    elif nhalo > 0:
        print('found halos in:', label)
        IsAvailable = True
        # halo_index = int(halo_prog_index[T])
        # halo_index = int(0)
        halo_x = Halo_Coor[label][halo_index, 0]
        halo_y = Halo_Coor[label][halo_index, 1]
        halo_z = Halo_Coor[label][halo_index, 2]
        halo_r = Halo_R_200[label][halo_index]
        halo_pos = [halo_x, halo_y, halo_z]
        halo_nsub = Halo_NSub[label][halo_index]
        if halo_nsub == 0:
            print ('no subhalo in:', label)
            continue
        print('Number of subhalos:', halo_nsub)
        
        halo_np = Halo_ParN[label][halo_index]
        halo_off = Halo_Offset[label][halo_index]
        # halo_off = int(0)
        
        par_x_inhalo = Par_x[int(halo_off):int(halo_off+halo_np)]
        par_y_inhalo = Par_y[int(halo_off):int(halo_off+halo_np)]
        par_z_inhalo = Par_z[int(halo_off):int(halo_off+halo_np)]
        par_id_inhalo = Par_id[int(halo_off):int(halo_off+halo_np)]
        
        par_selection = func_BoxSelect(Par_x, Par_y, Par_z,
                                    np.min(par_x_inhalo), np.max(par_x_inhalo), 
                                    np.min(par_y_inhalo), np.max(par_y_inhalo), 
                                    np.min(par_z_inhalo), np.max(par_z_inhalo))
        
        par_x_inbox = Par_x[par_selection]
        par_y_inbox = Par_y[par_selection]

        N_panel = np.min([max_npeak, halo_nsub, 10])
        # N_try = np.min([max_npeak, halo_nsub, 20])
        nrows = 4
        ncols = int(N_panel + 1)
        fig, axs = plt.subplots(figsize = (ncols*2+4, nrows*2+4),  #sharex = 'row', sharey = 'row',
                            nrows = nrows, ncols = ncols, 
                            dpi = 200)

        # dens_map, xedges, yedges = np.histogram2d(par_x_inbox, par_y_inbox, bins=[NbinH, NbinH]) 
        # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        # axs[1, 0].imshow(np.transpose(np.log10(dens_map+dens_factorH)), 
        #                         extent=extent, origin='lower',
        #                         cmap=cm.twilight_shifted, interpolation="spline16")
        axs[1, 0].scatter(halo_x, halo_y, alpha = 1,  
                        marker = '+', s = 100, c = 'black')
        
        par_inhalo_selection = index_parent.get_indexer(par_id_inhalo)
        par_inhalo_selection = par_inhalo_selection[par_inhalo_selection >= 0]
        par_x_inhalo0 = Par_x0[par_inhalo_selection]
        par_y_inhalo0 = Par_y0[par_inhalo_selection]
        par_z_inhalo0 = Par_z0[par_inhalo_selection]
        print('{} particles, found {} in Lagrangian'.format(len(par_id_inhalo), len(par_x_inhalo0)))
        par_linrho_inhalo = Par_linrho[par_inhalo_selection]
        CoorAllHaloPar0[label] = [par_x_inhalo0, par_y_inhalo0, par_z_inhalo0]
        LinRhoAllHaloPar[label] = par_linrho_inhalo

        func_scatter(axs[0, 0], par_x_inhalo0, par_y_inhalo0, 
                    Nrand, 1, 0.4, 'grey')
        
        mean_linrho_halo = np.mean(par_linrho_inhalo)
        std_linrho_halo = np.std(par_linrho_inhalo)
        print(mean_linrho_halo, std_linrho_halo)
        # lowlim_linrho = mean_linrho_halo + std_linrho_halo
        collapse_lim = 1.686/(Redshift[label] + 1)*(Redshift[label] + 1)

        par_linrho_inhalo_temp = par_linrho_inhalo
        normal = colors.Normalize(vmin = 1, vmax = N_panel)
        peak_ids = []

        # spacing, filter_size_final, filter_size, npeak_filter, grid, id_grid, peak_coords,\
        # x_grid_coord, y_grid_coord, z_grid_coord = find_min_peak(par_x_inhalo0, par_y_inhalo0, par_z_inhalo0, 
        #                             par_linrho_inhalo, par_id_inhalo)
        # npeaks = peak_coords.shape[0]
        # Peaks_Properties[label] = [npeaks, filter_size, npeak_filter]
        # grid, id_grid, x_grid_coord, y_grid_coord, z_grid_coord, sizes = \
        #                             set_grid(par_x_inhalo0, par_y_inhalo0, 
        #                                      par_z_inhalo0, par_linrho_inhalo, par_id_inhalo)
        try:
            npeaks = int(np.max(NPeak_Test[label][1]))
            grid, id_grid, x_grid_coord, y_grid_coord, z_grid_coord, sizes = Grid_Sizes[label]
            final_list = Peak_List[label]
            print('Found {} peaks in a {} list'.format(npeaks, len(final_list[0, :])))
        except:
            print('No data found!', label)
            continue

        # if np.sum(sizes != grid_size) == 0:
        #     print('Found {} peaks with filter size = {} grid size'.format(npeaks, grid_size))
        # else:
        #     print('Error: grid not the same!')
        # rfs_ngrid = rfs/spacing
        # print('Minimum wavelength is {} grid length'.format(rfs_ngrid))

        for i in range(npeaks):

            label_peak = '{}P{}'.format(label, i)
            print('looking at {}th peak'.format(i))

            peak_maxlinrho_temp = grid[int(final_list[0, i]), int(final_list[1, i]), int(final_list[2, i])]

            x_maxlinrho_temp = x_grid_coord[int(final_list[0, i])]
            y_maxlinrho_temp = y_grid_coord[int(final_list[1, i])]
            z_maxlinrho_temp = z_grid_coord[int(final_list[2, i])]
            id_maxlinrho_temp = int(id_grid[int(final_list[0, i]), int(final_list[1, i]), int(final_list[2, i])])
            
            index_maxlinrho_temp = Par_id == id_maxlinrho_temp
            print(id_maxlinrho_temp)
            
            if (np.sum(index_maxlinrho_temp) < 1):
                print('Error in finding the corresponding peak at z', label)
                continue
            if (np.sum(index_maxlinrho_temp) > 1):
                print('More than one particle sharing ID?', label)
                print(np.sum(index_maxlinrho_temp))
                print(Par_id[index_maxlinrho_temp])

            # exit()

            peak_ids.append(id_maxlinrho_temp)
            if i < N_panel:
                axs[0, 0].scatter(x_maxlinrho_temp, y_maxlinrho_temp, alpha = 0.8,  
                                marker = 'x', s = 100, c = cm.rainbow(normal(i)))
                print(x_maxlinrho_temp, y_maxlinrho_temp)
                
            par_x_maxlinrho_temp = Par_x[index_maxlinrho_temp]
            par_y_maxlinrho_temp = Par_y[index_maxlinrho_temp]
            par_z_maxlinrho_temp = Par_z[index_maxlinrho_temp]
            if i < N_panel:
                axs[1, 0].scatter(par_x_maxlinrho_temp, par_y_maxlinrho_temp, alpha = 0.8, 
                                marker = 'x', s = 100, c = cm.rainbow(normal(i)))

            centre_maxlinrho_temp = [x_maxlinrho_temp, y_maxlinrho_temp, z_maxlinrho_temp]
            peak_maxlinrho_temp = [par_x_maxlinrho_temp, par_y_maxlinrho_temp, par_z_maxlinrho_temp]
            dist_to_maxlinrho_temp = distance(par_x_inhalo0, par_y_inhalo0, par_z_inhalo0, 
                                            centre_maxlinrho_temp)
            near_maxlinrho_selection_temp = dist_to_maxlinrho_temp < rfs
            x_near_maxlinrho_temp = par_x_inhalo0[near_maxlinrho_selection_temp]
            y_near_maxlinrho_temp = par_y_inhalo0[near_maxlinrho_selection_temp]
            z_near_maxlinrho_temp = par_z_inhalo0[near_maxlinrho_selection_temp]
            id_near_maxlinrho_temp = par_id_inhalo[near_maxlinrho_selection_temp]

            # par_x_inhalo_peak = par_x_inhalo[near_maxlinrho_selection_temp]
            # par_y_inhalo_peak = par_x_inhalo[near_maxlinrho_selection_temp]
            # par_z_inhalo_peak = par_x_inhalo[near_maxlinrho_selection_temp]

            x_near_maxlinrho_temp_c = x_near_maxlinrho_temp - x_maxlinrho_temp
            y_near_maxlinrho_temp_c = y_near_maxlinrho_temp - y_maxlinrho_temp
            z_near_maxlinrho_temp_c = z_near_maxlinrho_temp - z_maxlinrho_temp
            
            x_axis = dist_to_maxlinrho_temp[near_maxlinrho_selection_temp]
            y_axis = par_linrho_inhalo[near_maxlinrho_selection_temp]

            dist_to_maxlinrho_all = distance(Par_x0, Par_y0, Par_z0, 
                                            centre_maxlinrho_temp)
            near_maxlinrho_all_selection = dist_to_maxlinrho_all < rfs
            linrho_all = Par_linrho[near_maxlinrho_all_selection]    
            dist_all = dist_to_maxlinrho_all[near_maxlinrho_all_selection]   
            if i < N_panel:                       
                func_scatter(axs[1, i+1], dist_all, linrho_all, Nrand, 1, 0.5, c = cm.rainbow(normal(i)))

                func_scatter(axs[0, i+1], x_axis, y_axis, Nrand, 1, 0.5, c = cm.rainbow(normal(i)))

            y_axis_selection = y_axis > 0
            # if np.sum(y_axis_selection) == 0:
            #     continue
            x_axis = x_axis[y_axis_selection]
            y_axis = y_axis[y_axis_selection]
            if len(y_axis) < 20:
                print('Not enough particles for fitting')
                continue
            x_near_maxlinrho_temp_c = x_near_maxlinrho_temp_c[y_axis_selection]
            y_near_maxlinrho_temp_c = y_near_maxlinrho_temp_c[y_axis_selection]
            z_near_maxlinrho_temp_c = z_near_maxlinrho_temp_c[y_axis_selection]
            if i < N_panel:
                axs[0, i+1].axvline(x = rfs, color = 'blue', 
                            linewidth = 2, linestyle = 'solid', alpha = 0.8)

            dist_axis = np.linspace(np.min(x_axis), np.max(x_axis), Nbin)
            quad_coef = [0, 0]
            # quad_coef[0], quad_coef[1] = func_Fitqualin(x_axis, y_axis)
            quad_coef[0], quad_coef[1] = np.polyfit(x_axis**2, y_axis, 1)
            print('{}ith peak, quad = '.format(quad_coef))
            quad_fit = np.polyval(quad_coef, dist_axis**2)
            peak_r = np.sqrt(-quad_coef[-1]/quad_coef[0])
            print('Fitted coefficient: ', quad_coef)
            if i < N_panel:
                # axs[0, i+1].plot(dist_axis, quad_fit, color = 'black', linewidth = 2, linestyle = 'dashed')
                axs[0, i+1].axvline(x = peak_r, color = 'red', 
                            linewidth = 2, linestyle = 'solid', alpha = 0.8)
            
            ell_fit_coef, pcov = curve_fit(func_FitQua_ellipsoid, 
                                        [x_near_maxlinrho_temp_c, 
                                        y_near_maxlinrho_temp_c, 
                                        z_near_maxlinrho_temp_c, 
                                        np.ones(len(x_near_maxlinrho_temp_c))*quad_coef[1]], 
                                        y_axis, 
                                        p0 = [quad_coef[0], quad_coef[0], quad_coef[0], 0, 0, 0])

            ell_fit_matrix = np.zeros((3, 3))
            ell_fit_matrix[0, 0] = ell_fit_coef[0]
            ell_fit_matrix[1, 1] = ell_fit_coef[1]
            ell_fit_matrix[2, 2] = ell_fit_coef[2]
            ell_fit_matrix[0, 1] = ell_fit_matrix[1, 0] = ell_fit_coef[5]/2
            ell_fit_matrix[2, 1] = ell_fit_matrix[1, 2] = ell_fit_coef[3]/2
            ell_fit_matrix[0, 2] = ell_fit_matrix[2, 0] = ell_fit_coef[4]/2

            fitted_r = np.sqrt(func_FitQua_ellipsoid([x_near_maxlinrho_temp_c, 
                                                    y_near_maxlinrho_temp_c, 
                                                    z_near_maxlinrho_temp_c, 0], 
                                                *ell_fit_coef)/quad_coef[0])
            
            if i < N_panel:
                func_scatter(axs[2, i+1], fitted_r, y_axis, Nrand, 0.5, 0.5, c = 'grey')
                func_scatter(axs[3, i+1], fitted_r, y_axis, Nrand, 0.5, 0.5, c = 'grey')

            index_bin_min, dist_bin_min, linrho_bin_min, \
            index_bin_max, dist_bin_max, linrho_bin_max = func_minmax(fitted_r, y_axis, Ndistbin)
            quad_coef_upper = [0, 0]
            quad_coef_upper[0], quad_coef_upper[1] = func_Fitqualin(dist_bin_max, linrho_bin_max)

            par_x_inhalo0_c = par_x_inhalo0 - x_maxlinrho_temp
            par_y_inhalo0_c = par_y_inhalo0 - y_maxlinrho_temp
            par_z_inhalo0_c = par_z_inhalo0 - z_maxlinrho_temp
            dist_inhalo0_ell = np.sqrt(func_FitQua_ellipsoid([par_x_inhalo0_c, par_y_inhalo0_c, par_z_inhalo0_c, 0], 
                                                *ell_fit_coef)/quad_coef[0])
            quad_fit_upper_halo = np.polyval(quad_coef_upper, dist_inhalo0_ell**2)
            semi_selection = par_linrho_inhalo < quad_fit_upper_halo
            
            print('semi fraction:', np.sum(semi_selection)/len(semi_selection), label_peak)

            par_x_semi = par_x_inhalo[semi_selection]
            par_y_semi = par_y_inhalo[semi_selection]
            par_z_semi = par_z_inhalo[semi_selection]
            Semi_Selection[label_peak] = semi_selection


            eigenvals, eigenvecs = np.linalg.eig(ell_fit_matrix)
            eigenvals = np.sort(eigenvals**2)
            lambda1 = eigenvals[2]
            lambda2 = eigenvals[1]
            lambda3 = eigenvals[0]
            ellipticity = (lambda1 - lambda3)/2/(lambda1+lambda2+lambda3)
            prolateness = (lambda1 + lambda3 - 2*lambda2)/2/(lambda1+lambda2+lambda3)

            
            colors_list = ['salmon', 'teal', 'black', 'blue', 'red']
            # divider_list = [50, 10, 5, 2, 1]
            divider_list = [2, 1]
            # print(fitted_r)
            log_coefs = []

            dist_axis = np.linspace(np.min(fitted_r), np.max(fitted_r), 100)
            for dv_index, divider in enumerate(divider_list):
                dist_near_peak_fit_index = np.argsort(fitted_r)[:int(len(fitted_r)/divider)]
                dist_near_peak_fit = fitted_r[dist_near_peak_fit_index]
                linrho_near_peak_fit = y_axis[dist_near_peak_fit_index]
                if len(dist_near_peak_fit) > 100:
                    if i < N_panel:
                        for temp_irow in [2, 3]:
                            # func_scatter(axs[temp_irow, i+1], fitted_r, y_axis, Nrand, 0.5, 0.5, c = 'grey')
                            axs[temp_irow, i+1].axvline(x = np.max(dist_near_peak_fit), 
                                            color = colors_list[dv_index], linestyle = 'dashed')
                            axs[temp_irow, i+1].set_ylim(np.min(y_axis), np.max(y_axis))

                    quad_coef_dv = [0, 0]
                    quad_coef_dv[0], quad_coef_dv[1] = func_Fitqualin(dist_near_peak_fit, 
                                                                    linrho_near_peak_fit)
                    quad_fit_dv = np.polyval(quad_coef_dv, dist_axis**2)
                    if i < N_panel:
                        axs[2, i+1].plot(dist_axis, quad_fit_dv, 
                                color = colors_list[dv_index], linewidth = 1, 
                                linestyle = 'dashed')

                    sub_y_axis = np.abs(linrho_near_peak_fit - np.max(linrho_near_peak_fit))
                    y_axis_log = np.log10(sub_y_axis)
                    x_axis_log = np.log10(dist_near_peak_fit)
                    logx_axis_selection = dist_near_peak_fit > 0
                    logy_axis_selection = sub_y_axis > 0
                    log_selection = logx_axis_selection*logy_axis_selection
                    if np.sum(log_selection) < 5:
                        continue
                    y_axis_log = y_axis_log[log_selection]
                    x_axis_log = x_axis_log[log_selection]

                    log_coef = np.polyfit(x_axis_log, y_axis_log, 1)
                    print(log_coef)
                    log_coefs.append(log_coef)
                    poly3_fit_dv = 10**np.polyval( log_coef, np.log10(dist_axis))
                    if i < N_panel:
                        axs[3, i+1].plot(dist_axis, -poly3_fit_dv + np.max(linrho_near_peak_fit), 
                                color = colors_list[dv_index], linewidth = 1, 
                                linestyle = 'dashed')

                    # poly4_coef, pcov = curve_fit(func_4PowerLaw, dist_near_peak_fit, 
                    #                     linrho_near_peak_fit - np.max(linrho_near_peak_fit), 
                    #                     p0 = [0, poly3_coef[0], quad_coef_dv[0]])
                    # poly4_fit_dv = func_4PowerLaw(dist_axis, *poly4_coef)
                    # axs[4, i].plot(dist_axis, poly4_fit_dv + np.max(linrho_near_peak_fit), 
                    #             color = colors_list[dv_index], linewidth = 1, 
                    #             linestyle = 'dashed')

                    # exp_coef, pcov = curve_fit(func_2Exp, dist_near_peak_fit, 
                    #                     linrho_near_peak_fit - np.max(linrho_near_peak_fit))
                    # exp2_fit_dv = func_2Exp(dist_axis, *exp_coef)
                    # axs[5, i].plot(dist_axis, exp2_fit_dv + np.max(linrho_near_peak_fit), 
                    #             color = colors_list[dv_index], linewidth = 1, 
                    #             linestyle = 'dashed')
                    
            # inertia = func_inertia(x_near_maxlinrho_temp[y_axis_selection], 
            #                     y_near_maxlinrho_temp[y_axis_selection],
            #                     z_near_maxlinrho_temp[y_axis_selection])
            # eigenval_par, eigenvec_par =  lina.eigh(inertia)
            # ca = np.max(eigenval_par)/np.min(eigenval_par)

            label_peak = '{}P{}'.format(label, i)
            print('save data', label_peak)

            Peaks_Properties[label_peak] = [centre_maxlinrho_temp, peak_maxlinrho_temp, 
                                            x_axis, y_axis, 
                                            quad_coef, log_coefs, ell_fit_coef, 
                                            peak_r, ellipticity, prolateness, 
                                            eigenvals, eigenvecs, ell_fit_matrix,
                                            # ca,eigenvals
                                            id_near_maxlinrho_temp[y_axis_selection],
                                            x_near_maxlinrho_temp[y_axis_selection], 
                                            y_near_maxlinrho_temp[y_axis_selection], 
                                            z_near_maxlinrho_temp[y_axis_selection]]
            # Peaks_Properties[label] = i
            # near_maxlinrho_selection_peak = dist_to_maxlinrho_temp < peak_r

            # par_linrho_nearpeak_temp = par_linrho_inhalo_temp[near_maxlinrho_selection_peak]
            # frac_volme_linrho0 = np.sum(par_linrho_nearpeak_temp <= 1e-5)/np.sum(near_maxlinrho_selection_peak)
            # if frac_volme_linrho0 < 



            # par_linrho_inhalo_temp[near_maxlinrho_selection_peak] = 0
            # mean_linrho_halo_temp = np.mean(par_linrho_inhalo_temp)

            # if mean_linrho_halo_temp < 0:
            #     print('break at {}th peak'.format(i))
            #     
            #     break
            

        # print('Number of peaks {} in {} try'.format(i, N_try))
        Peak_IDs[label] = peak_ids
        for irow in range(nrows):
            for icol in range(ncols):
                axs[irow, icol].tick_params(axis='both', which='both', direction='in', top=True, right=True, length=6)
        plt.tight_layout()
        path = '/cosma7/data/dp004/dc-wang4/VVV_plots/Halo_Peak/'+\
        'AllPeak_{}.png'.format(label)
        print('save image: ', path)
        fig.savefig(path, bbox_inches = 'tight', dpi=300)
        plt.close(fig)

AllPeak_data = []
AllPeak_data.append(Peaks_Properties)
AllPeak_data.append(CoorAllHaloPar0)
AllPeak_data.append(LinRhoAllHaloPar)
AllPeak_data.append(Peak_IDs)
AllPeak_data.append(Semi_Selection)

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Peak_Data/AllPeak_{}'.format(label)
print('save data: ', path)
np.save(path, AllPeak_data)

#print('Done all: '+str((time.clock()-t1)/60)+' min')


        









