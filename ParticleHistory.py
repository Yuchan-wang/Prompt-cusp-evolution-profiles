import numpy as np
# import h5py as h5
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
# import pandas as pd
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
print('T min: ', T_min)
print('T max: ', T_max)
print('kfs: ', kfs)

f = open('/cosma7/data/dp004/dc-wang4/IC_gen/apply_ps/ps/extended_planck_linear_powspec', 'r')

lines = f.readlines()
nl = len(lines)
ks = np.zeros(nl)
delta = np.zeros(nl)
for i in range(nl):
    line = lines[i].strip('[]\n').split('  ')
    ks[i] = 10**float(line[0])
    delta[i] = np.sqrt(10**float(line[1]))

# k_factor = (ks/kfs)**2
# D_factor = (1-2/3*k_factor)*np.exp(-k_factor)
# delta_sf = delta*D_factor
# logdelta = np.log10(delta_sf)
# ks_l = ks[logdelta > -2]
# logdelta_l = logdelta[logdelta > -2]
# logdelta2 = 2*logdelta_l
# max_delta = np.argmax(logdelta2)
# k_max = ks_l[max_delta]

# print('k_max: ', k_max)

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

Indexer_Halo = {}
Halo_Dist = {}
Dist_Axis = {}
Halo_Density = {}
Par_SubIndex = {}
Halo_Index = {}
for T_min_temp in range(0, T_min+1):
    # for T_max_temp in range(200, T_min_temp, -1):
    label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min_temp, T_min_temp)
    path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/AllSubRef_{}.npy'.format(label)
    try:
        AllHalo_data = np.load(path, allow_pickle = True)
        
        Indexer_Halo.update(AllHalo_data[0])
        # Halo_Dist.update(AllHalo_data[1])
        # Dist_Axis.update(AllHalo_data[2])
        # Halo_Density.update(AllHalo_data[3])
        Par_SubIndex.update(AllHalo_data[4])
        Halo_Index.update(AllHalo_data[11])
        print(path)
    except:
        pass

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
        # Halo_Mass[label] = File[13]*10**10/Hconst
        Halo_Coor[label] = File[14]/Hconst
        # Subhalo_Coor[label] = File[15]/Hconst
        Halo_R_200[label] = File[16]/Hconst
        Redshift[label] = File[17]
        Box_size[label] = File[18]/Hconst
        # Halo_Vel[label] = File[19]/Hconst
        # Subhalo_Spin[label] = File[20]

        # if T > 0:
        #     Subhalo_FirstProg[label] = Read_merger_tree(L, T, 16, Path)
    except:
        T_max = T-1
        break
    

print('finished reading', T_min, T_max)
if T_max < T_min:
    print('no file found')
    exit()

# label_tree = 'L{}HaloID{}k{}Nh{}'.format(L, HaloID, kfs, Nh)
# print('Loading merger tree information: ')
# SubProg_Index = {}
# prog_path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/Prog_Index_{}.npy'.format(label_tree)
# prog_data = np.load(prog_path, allow_pickle = True)
# SubProg_Index = prog_data.item()
# sub_prog_T = np.array(SubProg_Index[label_tree][0])
# sub_prog_index = np.array(SubProg_Index[label_tree][1])
# print(len(sub_prog_T), len(sub_prog_index))

halo_index = 0
rfactor = 12/7


Par_FirstIn_T = {}
Par_FirstSubIndex = {}
for T in range(T_min, T_max+1):  
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
    print(label)

    if len(Halo_Mass_200[label]) == 0:
        print('No Halo at ', label)
        continue

    try:
        halo_index = int(Halo_Index[label])
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
        # halo_index = int(0)
        halo_x = Halo_Coor[label][halo_index, 0]
        halo_y = Halo_Coor[label][halo_index, 1]
        halo_z = Halo_Coor[label][halo_index, 2]
        halo_r = Halo_R_200[label][halo_index]
        halo_pos = [halo_x, halo_y, halo_z]
        halo_nsub = Halo_NSub[label][halo_index]
        halo_off = Halo_Offset[label][halo_index]
        halo_np = Halo_ParN[label][halo_index]
        if halo_nsub == 0:
            print ('no subhalo in:', label)
            continue
        print('{} subhalos in halo'.format(halo_nsub))

        # par_selection = func_SphSelect(Par_x, Par_y, Par_z, halo_r,
        #                             halo_x, halo_y, halo_z)
    
        # par_x_inhalo = Par_x[par_selection]
        # par_y_inhalo = Par_y[par_selection]
        # par_z_inhalo = Par_z[par_selection]
        # par_id_inhalo = Par_id[par_selection]
        par_x_inhalo = Par_x[int(halo_off):int(halo_off+halo_np)]
        par_y_inhalo = Par_y[int(halo_off):int(halo_off+halo_np)]
        par_z_inhalo = Par_z[int(halo_off):int(halo_off+halo_np)]
        # par_m_inhalo = Par_m[int(halo_off):int(halo_off+halo_np)]
        par_id_inhalo = Par_id[int(halo_off):int(halo_off+halo_np)]

        par_inhalo_selection = np.zeros(len(par_id_inhalo)) + T
        par_insub_index = np.zeros(len(par_id_inhalo)) - 1

        for countT, T_ref in enumerate(range(int(T-1), 0, -1)):
            
            try:
                label_ref = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
                index_parent = Indexer_Halo[label_ref]
                [halo_nsub, par_sub_index] = Par_SubIndex[label_ref]
            except:
                print(T_ref)
                continue
            # 
            par_inhalo_selection_temp = index_parent.get_indexer(par_id_inhalo)
            isin = par_inhalo_selection_temp >= 0

            print(label_ref, np.sum(isin))
            par_inhalo_selection[isin] = T_ref

            par_insub_index[isin] = par_sub_index[par_inhalo_selection_temp][isin]

        Par_FirstIn_T[label] = par_inhalo_selection
        Par_FirstSubIndex[label] = par_insub_index

        print(label, len(Par_FirstSubIndex[label]))

AllParHistory_data = []
AllParHistory_data.append(Par_FirstIn_T)
AllParHistory_data.append(Par_FirstSubIndex)

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Particle_Data/AllParHistory_{}'.format(label)
print('save data: ', path)
np.save(path, AllParHistory_data)

print('Done all: '+str((time.clock()-t1)/60)+' min')
            