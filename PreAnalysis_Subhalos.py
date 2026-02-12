import numpy as np
import sys
# import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib import colors
# from scipy.ndimage.filters import gaussian_filter
# import matplotlib.patches as pat
# plt.style.use('dark_background')
# plt.style.use('default')
import warnings; warnings.simplefilter('ignore')
cmap=cm.viridis
# from scipy.optimize import curve_fit
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
# import os
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

print('Level: ', L)
print('HaloID: ', HaloID)
print('Nh: ', Nh)
print('Nf: ', Nh)
print('T min: ', T_min)
print('T max: ', T_max)
print('kfs: ', kfs)

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
SubHalo_GrN = {}

Redshift = {}
Box_size = {}
Par_LinearDensity = {}
Subhalo_FirstProg = {}

Peaks_Properties = {}
LinRhoAllHaloPar = {}
Peak_IDs = {}

print('finish initializing')

Nrand = int(1e3)
dens_factorH = 100
Nbin = 100
Ndistbin = 100
Nrhobin = 100
N_slice = 0
Omega_m = 0.307
Omega_L = 0.693

# Halo_PeakCoor = {}
# Halo_PeakID = {}
# Halo_Peak_Rs = {}
# Halo_Mr_Peak = {}
# Halo_PeakLin_bin = {}
# Halo_meanlinrho_Peak = {}
# Halo_EllParams = {}
# Halo_PeakEllR = {}
# Halo_ElliFit = {}
# Halo_PeakQua = {}
# Halo_PeakElli = {}
# Halo_PeakProl = {}
# Halo_PeakPar = {}
# Sub_Info = {}
# Par_Peak = {}

Indexer_Halo = {}
Halo_Dist = {}
Dist_Axis = {}
Halo_Density = {}
Par_SubIndex = {}
Halo_Centre = {}
Halo_ParCoor = {}


r_high_lim = halor
N_nfw = 51
Density_AllSub_Real = {}
Density_AllSub_Spurious = {}
        
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

    
    try:
        File = Read_file(L, T, 16, Path)
        label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
        print(Path, label)

        Par_Coor[label] = File[0]/Hconst
        Par_Mass[label] = File[1]*10**10/Hconst
        Par_ID[label] = File[2]
        Par_Vel[label] = File[3]/Hconst
        Halo_Mass_200[label] = File[4]*10**10/Hconst
        SubHalo_Mass[label] = File[5]*10**10/Hconst
        Halo_FirSub[label] = File[6]
        Halo_NSub[label] = File[7]
        SubHalo_GrN[label] = File[8]
        Halo_Offset[label] = File[9]
        Halo_ParN[label] = File[10]
        SubHalo_Offset[label] = File[11]
        SubHalo_ParN[label] = File[12]
        Halo_Mass[label] = File[13]*10**10/Hconst
        Halo_Coor[label] = File[14]/Hconst
        Subhalo_Coor[label] = File[15]/Hconst
        Halo_R_200[label] = File[16]/Hconst
        Redshift[label] = File[17]
        Box_size[label] = File[18]/Hconst
        Halo_Vel[label] = File[19]/Hconst
        Subhalo_Spin[label] = File[20]

        if T > 0:
            Subhalo_FirstProg[label] = Read_merger_tree(L, T, 16, Path)
    except:
        T_max = T-1
        break

print('finished reading', T_min, T_max)
if T_max < T_min:
    print('no file found')
    exit()

label_tree = 'L{}HaloID{}k{}Nh{}'.format(L, HaloID, kfs, Nh)
print('Loading merger tree information: ')
SubProg_Index = {}
prog_path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/Prog_Index_{}.npy'.format(label_tree)
prog_data = np.load(prog_path, allow_pickle = True)
SubProg_Index = prog_data.item()
sub_prog_T = np.array(SubProg_Index[label_tree][0])
sub_prog_index = np.array(SubProg_Index[label_tree][1])
print(len(sub_prog_T), len(sub_prog_index))

mean_density = np.sum(Par_Mass[label])/Box_size[label]**3


Subref_Info = []
Halo_Index = {}
for T_ref in range(T_min, T_max+1):  


    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
    z_current = Redshift[label]
    Ez = np.sqrt(Omega_m * (1 + z_current)**3 + Omega_L)
    rho_crit = mean_density * Ez**2 / Omega_m

    print(label)
    if len(Halo_Mass_200[label]) == 0:
        print('No Halo at ', label)
        continue

    T_index_selection = sub_prog_T == T_ref
    if np.sum(T_index_selection) == 0:
        print('No progenitor at ', label)
        continue

    sub_to_halo = SubHalo_GrN[label]
    halo_index = int(sub_to_halo[sub_prog_index[T_index_selection][0]])
    Halo_Index[label] = halo_index

    print('The progenitor of the halo is index {} at T = {}'.format(halo_index, T_ref))
    
    par_id_ref = Par_ID[label]

    halo_x = Halo_Coor[label][halo_index, 0]
    halo_y = Halo_Coor[label][halo_index, 1]
    halo_z = Halo_Coor[label][halo_index, 2]
    halo_pos = [halo_x, halo_y, halo_z]
    Halo_Centre[label] = halo_pos
    Par_x = Par_Coor[label][:, 0]
    Par_y = Par_Coor[label][:, 1]
    Par_z = Par_Coor[label][:, 2]
    Par_m = Par_Mass[label]
    halo_r = Halo_R_200[label][halo_index]
    halo_nsub = Halo_NSub[label][halo_index]
    halo_fsub = Halo_FirSub[label][halo_index]
    halo_off = Halo_Offset[label][halo_index]
    halo_np = Halo_ParN[label][halo_index]


    if halo_nsub == 0:
        print('No subhalo at ', label)
        continue

    # par_selection = func_SphSelect(Par_x, Par_y, Par_z, halo_r*2,
    #                             halo_x, halo_y, halo_z)
    
    # par_x_inhalo = Par_x[par_selection]
    # par_y_inhalo = Par_y[par_selection]
    # par_z_inhalo = Par_z[par_selection]
    # par_m_inhalo = Par_m[par_selection]
    # par_id_inhalo = par_id_ref[par_selection]
    par_x_inhalo = Par_x[int(halo_off):int(halo_off+halo_np)]
    par_y_inhalo = Par_y[int(halo_off):int(halo_off+halo_np)]
    par_z_inhalo = Par_z[int(halo_off):int(halo_off+halo_np)]
    par_m_inhalo = Par_m[int(halo_off):int(halo_off+halo_np)]
    par_id_inhalo = par_id_ref[int(halo_off):int(halo_off+halo_np)]
    Indexer_Halo[label] = pd.Index(par_id_inhalo)
    print(int(halo_off), int(halo_off+halo_np), halo_np)

    dist_halo = distance(par_x_inhalo, par_y_inhalo, par_z_inhalo, 
                        [halo_x, halo_y, halo_z])
    Halo_Dist[label] = dist_halo
    Halo_ParCoor[label] = [par_x_inhalo, par_y_inhalo, par_z_inhalo]

    r_conv = power_radius(par_x_inhalo, par_y_inhalo, par_z_inhalo, 
                        par_m_inhalo, halo_pos, rho_crit)
    dist_min = max(min(dist_halo), r_conv)
    # dist_max = max(dist_halo)
    dist_max = halo_r

    if dist_max > dist_min:

        dist_bin = 10**np.linspace(np.log10(dist_min), np.log10(dist_max), N_nfw)
        d_dist = dist_bin[1:] - dist_bin[:-1]
        dist_bin_mid = (dist_bin[1:] + dist_bin[:-1])/2
        # Dist_Axis[label] = [dist_bin, d_dist, dist_bin_mid]

        halo_count, subbin_edges = np.histogram(dist_halo, bins = dist_bin)
        density_halo = halo_count/d_dist/dist_bin_mid**2/4/np.pi*par_m_inhalo[0]
        # Halo_Density[label] = [dist_bin_mid, density_halo]

        r_low_ratio = dist_min / halo_r
        r_high_ratio = 1.0

        thetas, phis, dist_bin_mid_med, halo_dens_med, halo_dens_err, halo_dens_all = \
        func_HaloMedDensity(par_x_inhalo, par_y_inhalo, par_z_inhalo, 
                            par_m_inhalo, halo_r, halo_pos, 
                            r_low_ratio, r_high_ratio, N_nfw, 10)
        
        dist_bin_edges = 10**np.linspace(np.log10(dist_min), np.log10(dist_max), N_nfw+1)
        d_dist = dist_bin_edges[1:] - dist_bin_edges[:-1]
        dist_bin_mid_rec = (dist_bin_edges[1:] + dist_bin_edges[:-1])/2

        Dist_Axis[label] = [dist_bin_edges, d_dist, dist_bin_mid_rec]
        Halo_Density[label] = [dist_bin_mid_med, halo_dens_med]

    par_sub_index = np.zeros(len(par_x_inhalo)) - 1
    for isub in range(halo_nsub):
        
        print('looking at {}th subhalo'.format(isub))
        par_insub_off = SubHalo_Offset[label][int(isub)]
        par_insub_np = SubHalo_ParN[label][int(isub)]

        par_off_low = int(par_insub_off - halo_off)
        par_off_high = int(par_insub_off - halo_off + par_insub_np)
        print(int(par_insub_off), int(par_insub_off+par_insub_np), par_insub_np, 
            par_off_low, par_off_high)

        par_sub_index[par_off_low:par_off_high] = isub

    Par_SubIndex[label] = [halo_nsub, par_sub_index]

AllSub_data = []
AllSub_data.append(Indexer_Halo)
AllSub_data.append(Halo_Dist)
AllSub_data.append(Dist_Axis)
AllSub_data.append(Halo_Density)
AllSub_data.append(Par_SubIndex)
AllSub_data.append(Subhalo_FirstProg)
AllSub_data.append(SubHalo_Mass)
AllSub_data.append(SubHalo_GrN)
AllSub_data.append(Halo_Mass)
AllSub_data.append(Halo_Centre)
AllSub_data.append(Halo_ParCoor)
AllSub_data.append(Halo_Index)



label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/AllSubRef_{}'.format(label)
print('save data: ', path)
np.save(path, AllSub_data)

print('Done all: '+str((time.clock()-t1)/60)+' min')
