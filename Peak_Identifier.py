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


t1=time.clock()

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
redshift0 = Redshift[label0]
        
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

# Halo_Index = {}
# label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_min)
# path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/AllSubRef_{}.npy'.format(label)
# try:
#     AllHalo_data = np.load(path, allow_pickle = True)
#     Halo_Index.update(AllHalo_data[11])
#     print(path)
# except:
#     pass

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
        print('Halo Index:', halo_index)
    except:
        print('No progenitor at ', label)
        continue

    Par_x = Par_Coor[label][:, 0]
    Par_y = Par_Coor[label][:, 1]
    Par_z = Par_Coor[label][:, 2]
    Par_id = Par_ID[label]

    nhalo = len(Halo_Mass_200[label])
    if nhalo == 0:
        T_min = T+1
        print('no halo in: ', label)
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
        
        par_inhalo_selection = index_parent.get_indexer(par_id_inhalo)
        par_inhalo_selection = par_inhalo_selection[par_inhalo_selection >= 0]

        par_x_inhalo0 = Par_x0[par_inhalo_selection]
        par_y_inhalo0 = Par_y0[par_inhalo_selection]
        par_z_inhalo0 = Par_z0[par_inhalo_selection]
        par_linrho_inhalo = Par_linrho[par_inhalo_selection]

        #set up a grid
        boxl_x = np.max(par_x_inhalo0) - np.min(par_x_inhalo0)
        boxl_y = np.max(par_y_inhalo0) - np.min(par_y_inhalo0)
        boxl_z = np.max(par_z_inhalo0) - np.min(par_z_inhalo0)

        n_high_mass = 100

        # Nfactor = 1

        Nfactor_x = Nfactor_y = Nfactor_z = 1
        Nfactor_list = np.array([Nfactor_x, Nfactor_y, Nfactor_z], dtype = 'float')
        Dfactor = 0.1
        Grid_iterations = 0
        while n_high_mass > 0:
            par_selection0 = func_BoxSelect(Par_x0, Par_y0, Par_z0,
                                        (np.max(par_x_inhalo0) + np.min(par_x_inhalo0))/2 - boxl_x/Nfactor_list[0]/2, 
                                        (np.max(par_x_inhalo0) + np.min(par_x_inhalo0))/2 + boxl_x/Nfactor_list[0]/2, 
                                        (np.max(par_y_inhalo0) + np.min(par_y_inhalo0))/2 - boxl_y/Nfactor_list[1]/2, 
                                        (np.max(par_y_inhalo0) + np.min(par_y_inhalo0))/2 + boxl_y/Nfactor_list[1]/2, 
                                        (np.max(par_z_inhalo0) + np.min(par_z_inhalo0))/2 - boxl_z/Nfactor_list[2]/2, 
                                        (np.max(par_z_inhalo0) + np.min(par_z_inhalo0))/2 + boxl_z/Nfactor_list[2]/2, )
        
            par_x_inhalo0 = Par_x0[par_selection0]
            par_y_inhalo0 = Par_y0[par_selection0]
            par_z_inhalo0 = Par_z0[par_selection0]
            par_linrho_inhalo = Par_linrho[par_selection0]
            par_id_inhalo0 = Par_id[par_selection0]
        
            n_high_mass = np.sum(Par_m0[par_selection0] > Par_m0[0])

            # Nfactor += 0.1
            boxl_list = np.array([boxl_x/Nfactor_list[0], boxl_y/Nfactor_list[1], boxl_z/Nfactor_list[2]])
            box_max_arg = np.argmax(boxl_list)
            Nfactor_list[box_max_arg] = Nfactor_list[box_max_arg] + Dfactor
            Grid_iterations += 1
            print('Region found with high mass particles!', Nfactor_list, Grid_iterations)
            
            if np.sum(par_selection0) < 25**3:
                print('No region can be found without high mass particles!', Nfactor_list, Grid_iterations)
                sys.exit()


        par_x0_sorted = np.sort(np.unique(par_x_inhalo0))
        spacing = np.max(par_x0_sorted[1:] - par_x0_sorted[:-1])  

        grid, id_grid, x_grid_coord, y_grid_coord, z_grid_coord, sizes = \
            set_grid(par_x_inhalo0, par_y_inhalo0, 
                     par_z_inhalo0, par_linrho_inhalo, par_id_inhalo0)
        print('Grid setup: ', sizes)
        print('Iterations', Grid_iterations)
        n_grids = sizes[0]*sizes[1]*sizes[2]

        N_max_search = 10000
        r_search = 1

        collapse_lim = 1.686/(redshift0 + 1)*(Redshift[label] + 1)


        Npoints_test_list = range(int(n_grids/200), int(n_grids/40), int(n_grids/200))
        final_list_n = np.zeros(len(Npoints_test_list))

        for iN, N_points in enumerate(Npoints_test_list):
            final_list= find_max_grid(grid, N_points, r_search, collapse_lim, N_max_search)

            final_list = final_list[:, ~np.isnan(final_list[0])]
            final_list = np.unique(final_list, axis=1)
            
            final_list_n[iN] = len(final_list[0, :])

        NPeak_Test[label] = [Npoints_test_list, final_list_n]

        if len(final_list[0, :]) >= np.max(final_list_n):
            Peak_List[label] = final_list
        Grid_Sizes[label] = [grid, id_grid, x_grid_coord, y_grid_coord, z_grid_coord, sizes]

        

AllPeak_data = []
AllPeak_data.append(NPeak_Test)
AllPeak_data.append(Peak_List)
AllPeak_data.append(Grid_Sizes)

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Peak_Data/Identifier_{}'.format(label)
print('save data: ', path)
np.save(path, AllPeak_data)

print('Done all: '+str((time.clock()-t1)/60)+' min')


        









