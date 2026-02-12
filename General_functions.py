import numpy as np
import h5py as h5
import sys
# import matplotlib.pyplot as plt
# from scipy.ndimage.filters import gaussian_filter
import warnings; warnings.simplefilter('ignore')
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
from collections import Counter
# from sklearn.decomposition import nearest_lparticle
# from sklearn.preprocessing import StandardScaler
import numpy.linalg as lina
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
sys.path.append('/cosma7/data/dp004/dc-wang4/Nexus/NEXUS_1.1/python/python_import/')
# import MMF
# from sklearn.neighbors import KernelDensity
# from matplotlib.ticker import FormatStrFormatter
# import pandas as pd
# from scipy.stats.stats import spearmanr  
# from scipy import signal
# from scipy import interpolate as interp
from scipy.signal import savgol_filter
# from scipy.stats.stats import pearsonr  
# import time
# from skimage.feature import peak_local_max

def Read_file_ini(path, ID_total, Linear_Density):

    for ID in range(ID_total):
        path_P = path + '.{}.hdf5'.format(ID)        
        f_Par = h5.File(path_P, 'r')

        if ID==0:
            redshift = f_Par['Header'].attrs['Redshift']
            Box_size = f_Par['Header'].attrs['BoxSize']
            Hubble = f_Par['Header'].attrs['BoxSize']
            N_par = f_Par['Header'].attrs['NumPart_Total']
            N_par_Type1 = f_Par['Header'].attrs['NumPart_Total'][0]
            
            Par_Coor = np.zeros((np.sum(N_par), 3), dtype=f_Par['PartType1/Coordinates'].dtype)
            Par_Mass = np.zeros(np.sum(N_par), dtype=f_Par['PartType1/Masses'].dtype)
            Par_Vel = np.zeros((np.sum(N_par), 3), dtype=f_Par['PartType1/Velocities'].dtype)
            
            Par_ID = np.zeros(np.sum(N_par), dtype=f_Par['PartType1/ParticleIDs'].dtype)
            if Linear_Density is True:
                Par_LDensity = np.zeros(np.sum(N_par), dtype=f_Par['PartType1/lin-overdensity'].dtype)

            N_par_TypeN_count = N_par[1]
            N_par_Type1_count = 0
            
        N_par_in = f_Par['Header'].attrs['NumPart_ThisFile']                      
           
        
        for Partyle_index in range(1, 6):
            
            low_TypeN = int(N_par_TypeN_count)
            low_Type1 = int(N_par_Type1_count)
            if N_par_in[Partyle_index] == 0:
                continue
                
            if Partyle_index == 1:
                low_parN = low_Type1
            else:
                low_parN = low_TypeN
                
            high_parN = int(low_parN + N_par_in[Partyle_index])
            Par_Coor[low_parN:high_parN, :] = np.array(f_Par['PartType{}/Coordinates'.format(Partyle_index)][:])
            Par_Mass[low_parN:high_parN] = np.array(f_Par['PartType{}/Masses'.format(Partyle_index)][:])
            Par_ID[low_parN:high_parN] = np.array(f_Par['PartType{}/ParticleIDs'.format(Partyle_index)][:])
            Par_Vel[low_parN:high_parN, :] = np.array(f_Par['PartType{}/Velocities'.format(Partyle_index)][:])
            
            if Linear_Density is True:
                Par_LDensity[low_parN:high_parN] = np.array(f_Par['PartType{}/lin-overdensity'.format(Partyle_index)][:])
                
            if Partyle_index == 1:
                N_par_Type1_count = high_parN
            else:
                N_par_TypeN_count = high_parN

        f_Par.close()
    
    if Linear_Density is True:
        return Par_Coor, Par_Mass, Par_ID, Par_Vel, Par_LDensity, redshift
    else:
        return Par_Coor, Par_Mass, Par_ID, Par_Vel, redshift

def Read_file(L, T, ID_total, Path):
#     Par_Coor = np.empty((0,3), dtype='<f4')
#     Par_Mass = np.empty((0), dtype='<f4')
#     Par_ID = np.empty((0), dtype='<u8')
#     Par_Vel = np.empty((0,3), dtype='<f4')
    
#     Par_Coor_other = np.empty((0,3), dtype='<f4')
#     Par_Mass_other = np.empty((0), dtype='<f4')
#     Par_ID_other = np.empty((0), dtype='<u8')
#     Par_Vel_other = np.empty((0,3), dtype='<f4')
    
    Halo_Mass_200 = np.empty((0), dtype='<f4')
    Halo_Coor = np.empty((0,3), dtype='<f4')
    Halo_Vel = np.empty((0,3), dtype='<f4')
    Halo_R_200 = np.empty((0), dtype='<f4')
    Halo_NSub = np.empty((0), dtype='<i4')
    Halo_Offset = np.empty((0), dtype='<i8')
    Halo_ParN = np.empty((0), dtype='<i8')
    Halo_FirSub = np.empty((0), dtype='<i8')
    Halo_Mass = np.empty((0), dtype='<f4')
    SubHalo_Mass = np.empty((0), dtype='<f4')    
    SubHalo_GrN = np.empty((0), dtype='<i8')    
    SubHalo_Offset = np.empty((0), dtype='<i8')
    SubHalo_ParN = np.empty((0), dtype='<i4')
    Subhalo_IDMostbound = np.empty((0), dtype='<u8')    
    Subhalo_Coor = np.empty((0,3), dtype='<f4')
    Subhalo_HalfmassRad = np.empty((0), dtype='<f4')
    Subhalo_Spin = np.empty((0,3), dtype='<f4')
    Subhalo_Vel = np.empty((0,3), dtype='<f4')
    Subhalo_VelDisp = np.empty((0), dtype='<f4')
    Subhalo_Vmax = np.empty((0), dtype='<f4')
    Subhalo_VmaxRad = np.empty((0), dtype='<f4')

    N_halo_count = 0
    N_subhalo_count = 0
    FirstFileSub = True
    for ID in range(ID_total):
        if T < 100:
            path_P = Path + 'snapdir_0{1:02}/snap_0{2:02}.{3}.hdf5'.format(L, T, T, ID)  
            path_H = Path + 'groups_0{1:02}/fof_subhalo_tab_0{1:02}.{3}.hdf5'.format(L, T, T, ID)
        else:
            path_P = Path + 'snapdir_{1:02}/snap_{2:02}.{3}.hdf5'.format(L, T, T, ID)  
            path_H = Path + 'groups_{1:02}/fof_subhalo_tab_{1:02}.{3}.hdf5'.format(L, T, T, ID)
        f_Par = h5.File(path_P, 'r')

        if ID==0:
            redshift = f_Par['Header'].attrs['Redshift']
            Box_size = f_Par['Header'].attrs['BoxSize']
            Hubble = f_Par['Header'].attrs['BoxSize']
            N_par = f_Par['Header'].attrs['NumPart_Total']
            N_par_Type1 = f_Par['Header'].attrs['NumPart_Total'][0]
            
            for partype in range(0, 6):
                try:
                    Par_Coor = np.zeros((np.sum(N_par), 3), dtype=f_Par['PartType{}/Coordinates'.format(partype)].dtype)
                    Par_Mass = np.zeros(np.sum(N_par), dtype=f_Par['PartType{}/Masses'.format(partype)].dtype)
                    Par_ID = np.zeros(np.sum(N_par), dtype=f_Par['PartType{}/ParticleIDs'.format(partype)].dtype)
                    Par_Vel = np.zeros((np.sum(N_par), 3), dtype=f_Par['PartType{}/Velocities'.format(partype)].dtype)
                    break
                except:
                    continue
            
            N_par_TypeN_count = N_par[1]
            N_par_Type1_count = 0
            
        N_par_in = f_Par['Header'].attrs['NumPart_ThisFile']                      
           
        
        for Partyle_index in range(1, 6):
            
            low_TypeN = int(N_par_TypeN_count)
            low_Type1 = int(N_par_Type1_count)
            if N_par_in[Partyle_index] == 0:
                continue
                
            if Partyle_index == 1:
                low_parN = low_Type1
            else:
                low_parN = low_TypeN
                
            high_parN = int(low_parN + N_par_in[Partyle_index])
            Par_Coor[low_parN:high_parN, :] = np.array(f_Par['PartType{}/Coordinates'.format(Partyle_index)][:])
            if L != 'VVV_L0':
                Par_Mass[low_parN:high_parN] = np.array(f_Par['PartType{}/Masses'.format(Partyle_index)][:])
            Par_ID[low_parN:high_parN] = np.array(f_Par['PartType{}/ParticleIDs'.format(Partyle_index)][:])
            Par_Vel[low_parN:high_parN, :] = np.array(f_Par['PartType{}/Velocities'.format(Partyle_index)][:])
                
            if Partyle_index == 1:
                N_par_Type1_count = high_parN
            else:
                N_par_TypeN_count = high_parN
        f_Par.close()
        
        
        
        f_H = h5.File(path_H, 'r')
        # print(path_H)
        if ID == 0:
            N_halo = f_H['Header'].attrs['Ngroups_Total']
            N_halo_in = f_H['Header'].attrs['Ngroups_ThisFile']
            if N_halo > 0:
                Group = f_H['Group']
                Halo_Mass_200 = np.zeros(N_halo, dtype=Group['Group_M_Crit200'].dtype)
                Halo_Mass = np.zeros(N_halo, dtype=Group['GroupMass'].dtype)
                Halo_Coor = np.zeros((N_halo,3), dtype=Group['GroupPos'].dtype)
                Halo_R_200 = np.zeros(N_halo, dtype=Group['Group_R_Crit200'].dtype)
                Halo_FirSub = np.zeros(N_halo,  dtype=Group['GroupFirstSub'].dtype)
                Halo_NSub = np.zeros(N_halo,  dtype=Group['GroupNsubs'].dtype)
                Halo_Offset = np.zeros(N_halo, dtype=Group['GroupOffsetType'].dtype)
                Halo_ParN = np.zeros(N_halo,  dtype=Group['GroupLen'].dtype)
                Halo_Vel = np.zeros((N_halo,3),  dtype=Group['GroupVel'].dtype)
                
        N_halo_in = f_H['Header'].attrs['Ngroups_ThisFile']
        if N_halo_in > 0:
            Group = f_H['Group']
            low_Halo = int(N_halo_count)
            high_Halo = int(low_Halo + N_halo_in)
            Halo_Mass_200[low_Halo:high_Halo] = np.array(Group['Group_M_Crit200'][:])
            Halo_Mass[low_Halo:high_Halo] = np.array(Group['GroupMass'][:])
            Halo_Coor[low_Halo:high_Halo, :] = np.array(Group['GroupPos'][:])
            Halo_Vel[low_Halo:high_Halo, :] = np.array(Group['GroupVel'][:])
            Halo_R_200[low_Halo:high_Halo] = np.array(Group['Group_R_Crit200'][:])
            Halo_FirSub[low_Halo:high_Halo] = np.array(Group['GroupFirstSub'][:])
            Halo_NSub[low_Halo:high_Halo] = np.array(Group['GroupNsubs'][:])
            Halo_Offset[low_Halo:high_Halo] = np.array(Group['GroupOffsetType'][:, 1])
            Halo_ParN[low_Halo:high_Halo] = np.array(Group['GroupLen'][:])
            N_halo_count = high_Halo
        
        N_subhalo_in = f_H['Header'].attrs['Nsubhalos_ThisFile']
        N_subhalo = f_H['Header'].attrs['Nsubhalos_Total']
        if FirstFileSub == True:     
            if N_subhalo_in > 0:
                Subhalo = f_H['Subhalo']
                SubHalo_Mass = np.zeros(N_subhalo, dtype=Subhalo['SubhaloMass'].dtype)
                SubHalo_GrN = np.zeros(N_subhalo, dtype=Subhalo['SubhaloGroupNr'].dtype)
                SubHalo_Offset = np.zeros(N_subhalo, dtype=Subhalo['SubhaloOffsetType'].dtype)
                SubHalo_ParN = np.zeros(N_subhalo,dtype=Subhalo['SubhaloLenType'].dtype)
                Subhalo_Coor = np.zeros((N_subhalo,3), dtype=Subhalo['SubhaloPos'].dtype)
                Subhalo_Spin = np.zeros((N_subhalo,3), dtype=Subhalo['SubhaloSpin'].dtype)
                Subhalo_Vel = np.zeros((N_subhalo,3),dtype=Subhalo['SubhaloVel'].dtype)
                FirstFileSub = False
                
        N_subhalo_in = f_H['Header'].attrs['Nsubhalos_ThisFile']
        if N_subhalo_in>0:
            Subhalo = f_H['Subhalo']
            low_Subhalo = int(N_subhalo_count)
            high_Subhalo = int(low_Subhalo + N_subhalo_in)
            SubHalo_Mass[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloMass'][:])
            SubHalo_GrN[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloGroupNr'][:])
            SubHalo_Offset[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloOffsetType'][:, 1])
            SubHalo_ParN[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloLenType'][:, 1])
            Subhalo_Spin[low_Subhalo:high_Subhalo, :] = np.array(Subhalo['SubhaloSpin'][:])
            Subhalo_Coor[low_Subhalo:high_Subhalo, :] = np.array(Subhalo['SubhaloPos'][:])
            Subhalo_Vel[low_Subhalo:high_Subhalo, :] = np.array(Subhalo['SubhaloVel'][:])
            N_subhalo_count = high_Subhalo
            
        f_H.close()

    return Par_Coor, Par_Mass, Par_ID, Par_Vel, \
           Halo_Mass_200, SubHalo_Mass, Halo_FirSub, Halo_NSub, \
           SubHalo_GrN, Halo_Offset, Halo_ParN, SubHalo_Offset, SubHalo_ParN, \
           Halo_Mass, Halo_Coor, Subhalo_Coor, \
           Halo_R_200, redshift, Box_size, Halo_Vel, Subhalo_Spin

def Read_Group(L, T, ID_total, Path):
    Halo_Mass_200 = np.empty((0), dtype='<f4')
    Halo_Coor = np.empty((0,3), dtype='<f4')
    Halo_Vel = np.empty((0,3), dtype='<f4')
    Halo_R_200 = np.empty((0), dtype='<f4')
    Halo_NSub = np.empty((0), dtype='<i4')
    Halo_Offset = np.empty((0), dtype='<i8')
    Halo_ParN = np.empty((0), dtype='<i8')
    Halo_FirSub = np.empty((0), dtype='<i8')
    Halo_Mass = np.empty((0), dtype='<f4')
    SubHalo_Mass = np.empty((0), dtype='<f4')    
    SubHalo_GrN = np.empty((0), dtype='<i8')    
    SubHalo_Offset = np.empty((0), dtype='<i8')
    SubHalo_ParN = np.empty((0), dtype='<i4')
    Subhalo_IDMostbound = np.empty((0), dtype='<u8')    
    Subhalo_Coor = np.empty((0,3), dtype='<f4')
    Subhalo_HalfmassRad = np.empty((0), dtype='<f4')
    Subhalo_Spin = np.empty((0,3), dtype='<f4')
    Subhalo_Vel = np.empty((0,3), dtype='<f4')
    Subhalo_VelDisp = np.empty((0), dtype='<f4')
    Subhalo_Vmax = np.empty((0), dtype='<f4')
    Subhalo_VmaxRad = np.empty((0), dtype='<f4')

    N_halo_count = 0
    N_subhalo_count = 0
    FirstFileSub = True
    for ID in range(ID_total):
        if T < 100: 
            path_H = Path + 'groups_0{1:02}/fof_subhalo_tab_0{1:02}.{3}.hdf5'.format(L, T, T, ID)
        else:
            path_H = Path + 'groups_{1:02}/fof_subhalo_tab_{1:02}.{3}.hdf5'.format(L, T, T, ID)
    
        f_H = h5.File(path_H, 'r')
        # print(path_H)
        if ID == 0:
            N_halo = f_H['Header'].attrs['Ngroups_Total']
            N_halo_in = f_H['Header'].attrs['Ngroups_ThisFile']
            if N_halo > 0:
                Group = f_H['Group']
                Halo_Mass_200 = np.zeros(N_halo, dtype=Group['Group_M_Crit200'].dtype)
                Halo_Mass = np.zeros(N_halo, dtype=Group['GroupMass'].dtype)
                Halo_Coor = np.zeros((N_halo,3), dtype=Group['GroupPos'].dtype)
                Halo_R_200 = np.zeros(N_halo, dtype=Group['Group_R_Crit200'].dtype)
                Halo_FirSub = np.zeros(N_halo,  dtype=Group['GroupFirstSub'].dtype)
                Halo_NSub = np.zeros(N_halo,  dtype=Group['GroupNsubs'].dtype)
                Halo_Offset = np.zeros(N_halo, dtype=Group['GroupOffsetType'].dtype)
                Halo_ParN = np.zeros(N_halo,  dtype=Group['GroupLen'].dtype)
                Halo_Vel = np.zeros((N_halo,3),  dtype=Group['GroupVel'].dtype)
                
        N_halo_in = f_H['Header'].attrs['Ngroups_ThisFile']
        if N_halo_in > 0:
            Group = f_H['Group']
            low_Halo = int(N_halo_count)
            high_Halo = int(low_Halo + N_halo_in)
            Halo_Mass_200[low_Halo:high_Halo] = np.array(Group['Group_M_Crit200'][:])
            Halo_Mass[low_Halo:high_Halo] = np.array(Group['GroupMass'][:])
            Halo_Coor[low_Halo:high_Halo, :] = np.array(Group['GroupPos'][:])
            Halo_Vel[low_Halo:high_Halo, :] = np.array(Group['GroupVel'][:])
            Halo_R_200[low_Halo:high_Halo] = np.array(Group['Group_R_Crit200'][:])
            Halo_FirSub[low_Halo:high_Halo] = np.array(Group['GroupFirstSub'][:])
            Halo_NSub[low_Halo:high_Halo] = np.array(Group['GroupNsubs'][:])
            Halo_Offset[low_Halo:high_Halo] = np.array(Group['GroupOffsetType'][:, 1])
            Halo_ParN[low_Halo:high_Halo] = np.array(Group['GroupLen'][:])
            N_halo_count = high_Halo
        
        N_subhalo_in = f_H['Header'].attrs['Nsubhalos_ThisFile']
        N_subhalo = f_H['Header'].attrs['Nsubhalos_Total']
        if FirstFileSub == True:     
            if N_subhalo_in > 0:
                Subhalo = f_H['Subhalo']
                SubHalo_Mass = np.zeros(N_subhalo, dtype=Subhalo['SubhaloMass'].dtype)
                SubHalo_GrN = np.zeros(N_subhalo, dtype=Subhalo['SubhaloGroupNr'].dtype)
                SubHalo_Offset = np.zeros(N_subhalo, dtype=Subhalo['SubhaloOffsetType'].dtype)
                SubHalo_ParN = np.zeros(N_subhalo,dtype=Subhalo['SubhaloLenType'].dtype)
                Subhalo_Coor = np.zeros((N_subhalo,3), dtype=Subhalo['SubhaloPos'].dtype)
                Subhalo_Spin = np.zeros((N_subhalo,3), dtype=Subhalo['SubhaloSpin'].dtype)
                Subhalo_Vel = np.zeros((N_subhalo,3),dtype=Subhalo['SubhaloVel'].dtype)
                FirstFileSub = False
                
        N_subhalo_in = f_H['Header'].attrs['Nsubhalos_ThisFile']
        if N_subhalo_in>0:
            Subhalo = f_H['Subhalo']
            low_Subhalo = int(N_subhalo_count)
            high_Subhalo = int(low_Subhalo + N_subhalo_in)
            SubHalo_Mass[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloMass'][:])
            SubHalo_GrN[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloGroupNr'][:])
            SubHalo_Offset[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloOffsetType'][:, 1])
            SubHalo_ParN[low_Subhalo:high_Subhalo] = np.array(Subhalo['SubhaloLenType'][:, 1])
            Subhalo_Spin[low_Subhalo:high_Subhalo, :] = np.array(Subhalo['SubhaloSpin'][:])
            Subhalo_Coor[low_Subhalo:high_Subhalo, :] = np.array(Subhalo['SubhaloPos'][:])
            Subhalo_Vel[low_Subhalo:high_Subhalo, :] = np.array(Subhalo['SubhaloVel'][:])
            N_subhalo_count = high_Subhalo
            
        f_H.close()

    return Halo_Mass_200, SubHalo_Mass, Halo_FirSub, Halo_NSub, \
           SubHalo_GrN, Halo_Offset, Halo_ParN, SubHalo_Offset, SubHalo_ParN, \
           Halo_Mass, Halo_Coor, Subhalo_Coor, \
           Halo_R_200, Halo_Vel, Subhalo_Spin

def func_Halo_Mass(Mass, Nbin, IsLog):
    
    if len(Mass)==0 :
        return 0, 0
    elif len(Mass[Mass>0]) == 0:
        return 0, 0
    else:
        Mass_min = min(Mass[Mass>0])
        Mass_max = max(Mass[~np.isinf(Mass)])
        print(Mass_min, Mass_max)
        if IsLog == True:
            Mass_logBin = 10**np.linspace(np.log10(Mass_min - 1e-11), np.log10(Mass_max + 1e-11), Nbin)
        else:
            Mass_logBin = np.linspace(Mass_min, Mass_max, Nbin)

        Halo_N = np.zeros(Nbin-1)
        for i in range(0, Nbin-1):
            Halo_N[i] = np.sum(Mass>=Mass_logBin[i])

        Mass_axis = (Mass_logBin[1:] + Mass_logBin[:-1])/2
        
        return Mass_axis, Halo_N

    
def func_Selection(x, y, z, r, x_lowlim, x_higlim, y_lowlim, y_higlim, z_lowlim, z_higlim):
    z_s = np.logical_and(z>z_lowlim, z<z_higlim)
    x_s = np.logical_and(x>x_lowlim, x<x_higlim)
    y_s = np.logical_and(y>y_lowlim, y<y_higlim)
    s = z_s*x_s*y_s
    
    x_sed = x[s]
    y_sed = y[s]
    z_sed = z[s]
    r_sed = r[s]
    
    return x_sed, y_sed, z_sed, r_sed

def func_BoxSelect(x, y, z, x_lowlim, x_higlim, y_lowlim, y_higlim, z_lowlim, z_higlim):
    z_s = np.logical_and(z>z_lowlim, z<z_higlim)
    x_s = np.logical_and(x>x_lowlim, x<x_higlim)
    y_s = np.logical_and(y>y_lowlim, y<y_higlim)
    s = z_s*x_s*y_s
    
    return s

def func_SphSelect(xc, yc, zc, r, x, y, z):
    dis = (x-xc)**2 + (y-yc)**2 + (z-zc)**2
    Selection = dis < r**2
    return Selection

def unit_vector(vector):
#     """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
#     """ Returns the angle in radians between vectors 'v1' and 'v2'::

#             >>> angle_between((1, 0, 0), (0, 1, 0))
#             1.5707963267948966
#             >>> angle_between((1, 0, 0), (1, 0, 0))
#             0.0
#             >>> angle_between((1, 0, 0), (-1, 0, 0))
#             3.141592653589793
#     """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def power_radius(x, y, z, mass, centre, rho_crit):
    LHS = 0.6
    dist = distance(x, y, z, centre)
    indices = np.argsort(dist)
    mass, dist = mass[indices], dist[indices]
    Nr = np.array(range(np.size(mass))) + 1 # number of particles within r
    RHS = (np.sqrt(200) / 8) * np.sqrt((4 * np.pi * rho_crit) / (3 * mass[0]))
    RHS = RHS * np.sqrt(Nr[1:]) * dist[1:] ** 1.5 / np.log(Nr[1:])
    r_power = dist[np.argmin(np.abs(LHS-RHS))]

    return r_power

def PowerR_fromDistance(dist, mass_crit_ratio):
    LHS = 0.6
    indices = np.argsort(dist)
    dist = dist[indices]
    Nr = np.array(range(np.size(dist))) + 1 # number of particles within r
    # RHS = (np.sqrt(200) / 8) * np.sqrt((4 * np.pi * rho_crit) / (3 * mass[0]))
    RHS = (np.sqrt(200) / 8) * np.sqrt((4 * np.pi * mass_crit_ratio) / (3))
    RHS = RHS * np.sqrt(Nr[1:]) * dist[1:] ** 1.5 / np.log(Nr[1:])
    r_power = dist[np.argmin(np.abs(LHS-RHS))]
    return r_power

def distance(x, y, z, centre):
    dx = x - centre[0]
    dy = y - centre[1]
    dz = z - centre[2]
    return np.sqrt(dx**2 + dy**2 + dz**2)

def func_logNFW(r, rs, log10rho_s):
    return log10rho_s - np.log10(r/rs) - 2*np.log10(1 + r/rs)

def NFW_fitting(pars_in_x, pars_in_y, pars_in_z, par_in_m, halo_r,
                halo_pos, r_low_lim, r_high_lim, N):
    dist = np.sqrt((pars_in_x - halo_pos[0])**2 + \
                   (pars_in_y - halo_pos[1])**2 + \
                   (pars_in_z - halo_pos[2])**2 )
    
    dist_min = max(min(dist), halo_r*r_low_lim)
    dist_max = min(max(dist), halo_r*r_high_lim)
        
    if dist_min > dist_max:
        print('error')
    print(dist_min/halo_r, dist_max/halo_r)
    
    dist_bin = np.linspace(np.log10(dist_min), np.log10(dist_max), N)
    halo_dens_r = np.zeros(N-1)
    dist_bin_mid = 10**((dist_bin[1:] + dist_bin[:-1])/2)
    

    for j in range(N-1):
        d_dist_bin = 10**dist_bin[j+1] - 10**dist_bin[j]
        dist_select = np.logical_and(np.log10(dist) > dist_bin[j], 
                                     np.log10(dist) < dist_bin[j+1])
        
        mass_in_dist_bin = np.sum(par_in_m[dist_select])
        halo_dens_r[j] = mass_in_dist_bin/(dist_bin_mid[j])**2/d_dist_bin/4/np.pi
        
    halo_dens_r_non0 = halo_dens_r[halo_dens_r!=0]
    dist_bin_mid_non0 = dist_bin_mid[halo_dens_r!=0]
    
    try:
        popt, pcov = curve_fit(func_logNFW, dist_bin_mid_non0, np.log10(halo_dens_r_non0), )
        return dist_bin_mid, halo_dens_r, popt[0], popt[1], pcov
    except:
        return dist_bin_mid, halo_dens_r, 0, 0, 0
    
def Read_merger_tree(L, T, ID_total, Path):
    N_subhalo_count = 0
    subhalo_FirstProg = []
    for ID in range(ID_total):
        if T < 100:
            path_H = Path + 'groups_0{1:02}/subhalo_prog_0{1:02}.{3}.hdf5'.format(L, T, T, ID)
        else:
            path_H = Path + 'groups_{1:02}/subhalo_prog_{1:02}.{3}.hdf5'.format(L, T, T, ID)
        f_H = h5.File(path_H, 'r')
        try: 
            Subhalo = f_H['Subhalo']
        except: 
            continue
        if ID==0:
            N_subhalo = f_H['Header'].attrs['Nsubhalos_Total']
            subhalo_FirstProg = np.zeros(N_subhalo, dtype=Subhalo['FirstProgSubhaloNr'].dtype)
            
        N_subhalo_in = f_H['Header'].attrs['Nsubhalos_ThisFile']
        low_Subhalo = int(N_subhalo_count)
        high_Subhalo = int(low_Subhalo + N_subhalo_in)
        subhalo_FirstProg[low_Subhalo:high_Subhalo] = np.array(Subhalo['FirstProgSubhaloNr'][:])
        N_subhalo_count = high_Subhalo
        f_H.close()
    return subhalo_FirstProg
        
def func_HaloDensity(pars_in_x, pars_in_y, pars_in_z, par_in_m, halo_r,
                halo_pos, r_low_lim, r_high_lim, N):
    dist = np.sqrt((pars_in_x - halo_pos[0])**2 + \
                   (pars_in_y - halo_pos[1])**2 + \
                   (pars_in_z - halo_pos[2])**2 )
    
    dist_min = max(min(dist), halo_r*r_low_lim)
    dist_max = min(max(dist), halo_r*r_high_lim)
        
    if dist_min > dist_max:
        print('error')
    # print(dist_min/halo_r, dist_max/halo_r)
    
    dist_bin = np.linspace(np.log10(dist_min), np.log10(dist_max), N)
    halo_dens_r = np.zeros(N-1)
    dist_bin_mid = 10**((dist_bin[1:] + dist_bin[:-1])/2)
    
    for j in range(N-1):
        d_dist_bin = 10**dist_bin[j+1] - 10**dist_bin[j]
        dist_select = np.logical_and(np.log10(dist) > dist_bin[j], 
                                     np.log10(dist) < dist_bin[j+1])
        
        mass_in_dist_bin = np.sum(par_in_m[dist_select])
        halo_dens_r[j] = mass_in_dist_bin/(dist_bin_mid[j])**2/d_dist_bin/4/np.pi
        
    halo_dens_r_non0 = halo_dens_r[halo_dens_r!=0]
    dist_bin_mid_non0 = dist_bin_mid[halo_dens_r!=0]   
    
    return dist_bin_mid_non0, halo_dens_r_non0

def func_UniformThetaBin(Nbin):
    total = 2
    vol_bin = total/Nbin

    possibles = np.linspace(0, np.pi, Nbin*1000)
    solutions = np.zeros(Nbin+1)
    solutions[-1] = possibles[-1]
    ang = 0

    for i in range(Nbin):
        solutions[i] = ang
        
        mid_ang = (possibles + ang)/2
        d_ang = possibles - ang
        value = np.sin(mid_ang)*d_ang - vol_bin
        solution = possibles[np.argmin(np.abs(value))]
        
        ang = solution
    return solutions
    

def func_HaloMedDensity(pars_in_x, pars_in_y, pars_in_z, par_in_m, halo_r,
                halo_pos, r_low_lim, r_high_lim, Nbinr, NbinA):
    
    par_x_c = pars_in_x - halo_pos[0]
    par_y_c = pars_in_y - halo_pos[1]
    par_z_c = pars_in_z - halo_pos[2]
    
    XsqPlusYsq = par_x_c**2 + par_y_c**2
    dist = np.sqrt(XsqPlusYsq + par_z_c**2)
    thetas = np.arctan2(np.sqrt(XsqPlusYsq), par_z_c)
    phis = np.arctan2(par_y_c, par_x_c)
    
    dist_min = max(min(dist), halo_r*r_low_lim)
    dist_max = min(max(dist), halo_r*r_high_lim)
    dist_bin = np.linspace(np.log10(dist_min), np.log10(dist_max), Nbinr+1)
    dist_bin_mid = 10**((dist_bin[1:] + dist_bin[:-1])/2)
    d_dist_bin = 10**dist_bin[1:] - 10**dist_bin[:-1]
    
    Nbin_theta = NbinA
    cell_vol = 4*np.pi/NbinA**2
    # Nbin_phi = NbinA
    # theta_bins = np.linspace(0, np.pi, Nbin_theta+1)
    theta_bins = func_UniformThetaBin(Nbin_theta)
    # phi_bins = np.linspace(-np.pi, np.pi, Nbin_phi+1)
    d_theta = theta_bins[1:] - theta_bins[:-1]
    theta_bins_mid = (theta_bins[1:] + theta_bins[:-1])/2
    # d_phi = phi_bins[1:] - phi_bins[:-1]
    
    halo_dens_r = np.zeros(Nbinr)
    halo_dens_r_agbin = []
    halo_dens_r_error = np.zeros(Nbinr)

    for i in range(Nbinr):
        dist_selection = (np.log10(dist) >= dist_bin[i])*(np.log10(dist) < dist_bin[i+1])
        par_m_inrbin = par_in_m[dist_selection]
        thetas_inrbin = thetas[dist_selection]
        phis_inrbin = phis[dist_selection]
        # print('particles in r bin:', len(phis_inrbin))
        if np.sum(dist_selection) < 1:
            continue
        
        # halo_mass_inbin, xedges, yedges = np.histogram2d(thetas_inrbin, phis_inrbin, 
        #                                         bins=[theta_bins, phi_bins], 
        #                                         weights = par_m_inrbin, normed = False, 
        #                                         density = False) 
        # total_mass = np.sum(halo_mass_inbin)
        # total_mass_minus = total_mass - halo_mass_inbin
        # halo_dens_minus = total_mass_minus / ( (dist_bin_mid[i])**2 * d_dist_bin[i])
        # halo_dens_minus = halo_dens_minus / (4*np.pi - ((d_theta*np.sin(theta_bins_mid))[:, None])*d_phi[None, :])

        # halo_dens = halo_mass_inbin / ( (dist_bin_mid[i])**2 * d_dist_bin[i])
        # halo_dens = halo_dens / ((d_theta*np.sin(theta_bins_mid))[:, None]) / d_phi[None, :]

        # halo_dens = np.zeros((Nbin_theta, Nbin_phi))
        halo_dens = []
        for j in range(Nbin_theta):
            theta_selection = (thetas_inrbin < theta_bins[j+1])*(thetas_inrbin >= theta_bins[j])
            par_m_inbin = par_m_inrbin[theta_selection]
            thetas_inbin = thetas_inrbin[theta_selection]
            phis_inbin = phis_inrbin[theta_selection]
            # print('particles in theta bin:', np.sum(theta_selection))
            if np.sum(theta_selection) < 1:
                continue
            # d_theta = theta_bins[j+1] - theta_bins[j]
            # theta = (theta_bins[j+1] + theta_bins[j])/2
            Nbin_phi = int(np.rint(2*np.pi*(d_theta[j]*np.sin(theta_bins_mid[j]))/cell_vol))
            phi_bins = np.linspace(-np.pi, np.pi, Nbin_phi+1)
            d_phi = phi_bins[1:] - phi_bins[:-1]
            # print('Number of phi bins:', Nbin_phi)

            halo_mass_inbin, xedges = np.histogram(phis_inbin, bins=phi_bins, 
                                                weights = par_m_inbin, #normed = False, 
                                                density = False) 
            # dens_inbin = 
            total_mass = np.sum(halo_mass_inbin)
            vol_bin = (d_phi*d_theta[j]*np.sin(theta_bins_mid[j]))*\
                    ((dist_bin_mid[i])**2*d_dist_bin[i])

            halo_dens.extend(halo_mass_inbin/vol_bin)
            # total_mass_minus = total_mass - halo_mass_inbin     
            # halo_dens_minus = total_mass_minus / ( (dist_bin_mid[i])**2 * d_dist_bin[i])          
            # halo_dens.extend(halo_dens_minus / (d_theta[j]*np.sin(theta_bins_mid[j]))*(2*np.pi - d_phi))
        #     for k in range(Nbin_phi):
        #         phi_selection = (phis_inrbin < phi_bins[k+1])*(phis_inrbin >= phi_bins[k])
        #         d_phi = phi_bins[k+1] - phi_bins[k]

        #         selection = theta_selection*phi_selection

        #         mass_in_dist_bin = np.sum(selection)*par_m0
        #         halo_dens[j, k] = mass_in_dist_bin/(dist_bin_mid[i])**2/d_dist_bin[i]/4/np.pi/np.sin(theta)/d_theta/d_phi
        # halo_dens_r_agbin[i, :] = halo_dens
        halo_dens_r_agbin.extend(halo_dens)
        halo_dens = np.array(halo_dens)
        halo_dens = halo_dens.flatten()
        # print(halo_dens)
        halo_dens_r[i] = np.median(halo_dens)
        halo_dens_r_error[i] = np.std(halo_dens)
        if halo_dens_r[i] == 0:
            halo_dens_r[i] = np.mean(halo_dens)
        
    halo_dens_r_non0 = halo_dens_r[halo_dens_r!=0]
    dist_bin_mid_non0 = dist_bin_mid[halo_dens_r!=0]
    halo_dens_r_error_non0 = halo_dens_r_error[halo_dens_r!=0]
    # halo_dens_r_agbin = []
    return thetas, phis, dist_bin_mid_non0, halo_dens_r_non0, halo_dens_r_error_non0, halo_dens_r_agbin
    
def func_HaloCaustics(r, rho, window_length, polyorder):
    dr = np.log10(r[1]) - np.log10(r[0])
    logrho = np.log10(rho)
    
    ratio = savgol_filter(logrho, window_length=window_length, polyorder=polyorder, deriv=1, 
                          delta = dr, mode='nearest')
    
    return ratio

def func_FindPeak_Sub0(L, HaloID, kfs, Nh, T, T_min, T_max):

    for T_ref in range(T, T_min-1, -1):
        label_final = 0
        label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
        try:
            [halo_nsub, chosen_sub] = Halo_PeakSub[label]
            if chosen_sub == 0:
                return label
            elif chosen_sub > 0:
                if T_ref > T_min:
                    continue
                elif T_ref == T_min:
                    for T_ref in range(T_min, T_max+1):
                        try:
                            label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
                            [halo_nsub, chosen_sub] = Halo_PeakSub[label]
                            if chosen_sub == 0:
                                return label
                            else:
                                if T_ref < T_max:
                                    continue
                                elif T_ref == T_max:
                                    print('no peak in sub 0 avaialble')
                                    exit()
                        except:
                            if T_ref < T_max:
                                continue
                            elif T_ref == T_max:
                                print('no peak avaialble')
                                exit()
                
        except:
            if T_ref > T_min:
                continue
            elif T_ref == T_min:
                for T_ref in range(T_min, T_max+1):
                    try:
                        label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
                        [halo_nsub, chosen_sub] = Halo_PeakSub[label]
                        if chosen_sub == 0:
                            return label
                        else:
                            if T_ref < T_max:
                                continue
                            elif T_ref == T_max:
                                print('no peak in sub 0 avaialble')
                                exit()
                    except:
                        if T_ref < T_max:
                            continue
                        elif T_ref == T_max:
                            print('no peak avaialble')
                            exit()

def func_FitQua_ellipsoid(vec, a, b, c, f, g, h):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    B1 = vec[3]
    y = x**2*a + y**2*b + c*z**2 + f*2*y*z + g*2*x*z + h*2*x*y
    qua = y + B1
    
    return qua

def func_REllipsoid(vec, a, b, c, f, g, h):
    x = vec[0]
    y = vec[1]
    z = vec[2]    
    y = x**2*a + y**2*b + c*z**2 + f*2*y*z + g*2*x*z + h*2*x*y
    
    return y

def func_scatter(ax, x_axis, y_axis, Nrand, s, alpha, c):
    index_array = np.array(range(len(x_axis)))
    if Nrand > len(x_axis):
        ax.scatter(x_axis, y_axis, s = s, alpha = alpha, color = c)
    else:
        rand_index = np.random.choice(index_array, Nrand, replace=False)
        ax.scatter(x_axis[rand_index],
                             y_axis[rand_index], s = s, alpha = alpha, color = c)
    
    return 0

def func_ShellN(x, N):
    x = x[x>0]
    shell = np.zeros(len(x))
    # x_bins = 10**np.linspace(np.min(np.log10(x)), np.max(np.log10(x)), N + 1)
    x_bins_sort = np.argsort(x)
    x_bins = np.linspace(np.min(x_bins_sort), np.max(x_bins_sort), N + 1)
    for i in range(N):
        x_selection = (x_bins_sort > x_bins[i])*(x_bins_sort <= x_bins[i+1])
        shell[x_selection] = i
        
    return shell
def func_VelSep(vel_alongr, Nbin):
    hist, bin_edges = np.histogram(np.abs(vel_alongr), bins=Nbin)
    bin_mid = (bin_edges[1:] + bin_edges[:-1])/2

    gradient1 = np.gradient(hist, 1)

    gradient_pre = gradient1[:-1] < 0
    gradient_this = gradient1[1:] > 0
    
    selection = gradient_pre*gradient_this

    index_array = np.arange(1, len(bin_mid))
    try:
        index_ = np.min(index_array[selection])
    except:
        index_ = len(selection)

    return bin_mid, hist, bin_mid[index_]

def reduced_inertia(x_r, y_r, z_r, a, b, c, eigenvector):
    I = np.zeros((3, 3))
    transfer = np.transpose(eigenvector)

    N_par = len(x_r)
    coor_original = np.zeros((3, N_par))
    coor_original[0, :] = x_r
    coor_original[1, :] = y_r
    coor_original[2, :] = z_r
    coor_trans = np.matmul(transfer, coor_original)
    
    R2 = coor_trans[0, :]**2/a**2 + coor_trans[1, :]**2/b**2 + coor_trans[2, :]**2/c**2
    r_selection = R2!=0
    x_r = x_r[r_selection]
    y_r = y_r[r_selection]
    z_r = z_r[r_selection]
    R2 = R2[r_selection]
    
    I[0, 0] = np.sum( x_r**2 / R2)
    I[1, 1] = np.sum( y_r**2 / R2)
    I[2, 2] = np.sum( z_r**2 / R2)
    I[0, 1] = I[1, 0] = np.sum( x_r*y_r / R2)
    I[0, 2] = I[2, 0] = np.sum( x_r*z_r / R2)
    I[1, 2] = I[2 ,1] = np.sum( y_r*z_r / R2)
    return I

def inside_ellipse(x, y, z, a, b, c, eigenvector):
    transfer = np.transpose(eigenvector)
    N_par = len(x)
    coor_original = np.zeros((3, N_par))
    coor_original[0, :] = x
    coor_original[1, :] = y
    coor_original[2, :] = z
    coor_trans = np.matmul(transfer, coor_original)
    isin = (coor_trans[0, :]/a)**2 + (coor_trans[1, :]/b)**2 + (coor_trans[2, :]/c)**2
    return (isin <= 1), coor_trans

def halo_ellipse(halo_c, par_x, par_y, par_z, halo_r, r_ell):
    N_iteration = 30
    x = par_x - halo_c[0]
    y = par_y - halo_c[1]
    z = par_z - halo_c[2]
    volume = (halo_r*r_ell)**3
    
    a = b = c = halo_r
    eigenvector = np.zeros((3, 3))
    eigenvector[0, 0] = eigenvector[1, 1] = eigenvector[2, 2] = 1
    ss = qs = np.zeros(N_iteration)
    a_s = b_s = c_s = np.zeros(N_iteration)
    a_v = b_v = c_v = np.zeros((N_iteration, 3))
    
    for i in range(N_iteration):
        
        eli_selection, coor_transformed = inside_ellipse(x, y, z, a, b, c, eigenvector)
        
        x_eli = x[eli_selection]
        y_eli = y[eli_selection]
        z_eli = z[eli_selection]
        if np.sum(eli_selection)<50:
            break
        
        I = reduced_inertia(x_eli, y_eli, z_eli, a, b, c, eigenvector)
        eigenvalue, eigenvector =  lina.eigh(I)
        eigen_volume = np.sqrt(eigenvalue[0]*eigenvalue[1]*eigenvalue[2])
        eigenvalue = np.sqrt(eigenvalue)*(volume/eigen_volume)**(1/3)
        a = eigenvalue[0]
        b = eigenvalue[1]
        c = eigenvalue[2]
        ss[i] = a/c
        qs[i] = b/c
            
        if np.max( ( ( (ss[i] - ss[i-1])/ss[i] )**2, ((qs[i] - qs[i-1])/qs[i])**2) ) < 5*10**(-3) or i == N_iteration-1:
            direction = eigenvector
            # ca = np.min(eigenvalue)/np.max(eigenvalue)
            return direction, eigenvalue, eli_selection, x_eli, y_eli, z_eli, coor_transformed
        
        elif i==N_iteration-1:
            direction = eigenvector
            # ca = np.min(eigenvalue)/np.max(eigenvalue)
            return direction, eigenvalue, eli_selection, x_eli, y_eli, z_eli, coor_transformed

def func_spherical_angle(pars_in_x, pars_in_y, pars_in_z, dist, halo_pos):
    
    par_x_c = pars_in_x - halo_pos[0]
    par_y_c = pars_in_y - halo_pos[1]
    par_z_c = pars_in_z - halo_pos[2]
    
    XsqPlusYsq = par_x_c**2 + par_y_c**2
    thetas = np.arctan2(np.sqrt(XsqPlusYsq), par_z_c)
    phis = np.arctan2(par_y_c, par_x_c)

    return thetas, phis

def func_theoretical_rhor(r, R, A, p, mean_rho):
    return A*mean_rho*(r/R)**(p)

def func_theoretical_Amp(alpha, peak_height, peak_dev, density_coll, fec):
    return alpha*density_coll*peak_height**(9/4)*peak_dev**(-3/4)*fec**(-3/2)

def func_solve_fec(e, p):
    if p>e:
        return 1
    else:
        x = np.linspace(0, 2, 10000)
        y1 = x
        y2 = 1 + 0.47*(5*(e**2 - p*np.abs(p))*x**2)**(0.615)

        index = np.argwhere(np.diff(np.sign(y1 - y2))).flatten()

        return np.min(x[index])

def func_par_separation(x):
    keys = Counter(x).keys()
    unique = np.sort(np.array(list(keys)))
    separation = np.min(unique[1:] - unique[:-1])
    
    return separation

def func_minmax(x, y, N):
    y = y[x > 0]
    x = x[x > 0]
    # x_axis = 10**np.linspace(np.log10(np.min(x)), np.log10(np.max(x)), N+1)
    x_axis = np.linspace(np.min(x), np.max(x), N+1)
    x_axis_mid = (x_axis[1:] + x_axis[:-1])/2
    y_axis_min = np.zeros(N)
    y_axis_max = np.zeros(N)
    index_min = np.zeros(N)
    index_max = np.zeros(N)
    index_array = np.arange(0, len(x))
    IsEmpty = np.ones(N)
    for i in range(N):
        x_selection = (x < x_axis[i+1])*(x >= x_axis[i])
        y_selected = y[x_selection]
        index_selected = index_array[x_selection]
        if np.sum(x_selection) == 0:
            IsEmpty[i] = 0
            continue
        y_axis_min[i] = np.min(y_selected)
        index_min[i] = int(index_selected[int(np.argmin(y_selected))])
        y_axis_max[i] = np.max(y_selected)
        index_max[i] = int(index_selected[int(np.argmax(y_selected))])
        
    x_axis_mid_min = x_axis_mid[IsEmpty == 1]
    index_min = index_min[IsEmpty == 1]
    y_axis_min = y_axis_min[IsEmpty == 1]
    
    x_axis_mid_max = x_axis_mid[IsEmpty == 1]
    index_max = index_max[IsEmpty == 1]
    y_axis_max = y_axis_max[IsEmpty == 1]
    
    return index_min, x_axis_mid_min, y_axis_min, index_max, x_axis_mid_max, y_axis_max

def func_Fitqualin(x, y):
    y = y[x>0]
    x = x[x>0]
    maxy = np.max(y)
    y = y - maxy
    k = np.mean(y/x**2)
    
    return k, maxy

def func_3PowerLaw(x, a, b):
    y = a*x**3 + b*x**2
    return y

def func_4PowerLaw(x, a, b, c):
    y = a*x**4 + b*x**3 + c*x**2
    return y

def func_2Exp(x, a, b):
    y = a*np.exp(x*b)
    return y

def func_PowerLaw(x, a, b):
    y = a*x**b
    return y

def func_inertia(x, y, z):
    xc = np.mean(x)
    yc = np.mean(y)
    zc = np.mean(z)
    x_c = x - xc
    y_c = y - yc
    z_c = z - zc

    x2 = x_c**2 + y_c**2 + z_c**2

    I = np.zeros((3, 3))
    I[0, 0] = np.sum(y_c**2 + z_c**2)
    I[1, 1] = np.sum(z_c**2 + x_c**2)
    I[2, 2] = np.sum(y_c**2 + x_c**2)
    I[0, 1] = I[1, 0] = -np.sum(x_c*y_c)
    I[0, 2] = I[2, 0] = -np.sum(x_c*z_c)
    I[2, 1] = I[1, 2] = -np.sum(z_c*y_c)

    return I

def func_absdiff(y1, y2):
    return np.std(y1 - y2)

def func_fitR(r, rho, rfactor):
    rho_r = rho*r**rfactor
    amp = np.mean(rho_r)
    return amp
    
def func_fitNFW_rho(r, rho):
    popt, pcov = curve_fit(func_logNFW, r, np.log10(rho))
    return popt, pcov

def func_fitPower_rho(r, rho):
    popt, pcov = curve_fit(func_logPowerLaw, r, np.log10(rho))
    return popt, pcov

def func_logPowerLaw(r, a, b):
    return np.log10(a) + b*np.log10(r)

def Potential_Overlap(x1, y1, z1, x2, y2, z2, w1, w2, par_separation):
    
    m = len(x1)
    n = len(x2)
    
    R = 0
    for i in range(int(np.ceil(n/100))):
        low_index = i*100
        high_index = np.min([n, (i+1)*100])

        n_temp = high_index - low_index
        diff_mat = np.zeros((n_temp, m))
        x2_temp = x2[low_index:high_index]
        y2_temp = y2[low_index:high_index]
        z2_temp = z2[low_index:high_index]
        w2_temp = w2[low_index:high_index]
        
        w1_temp = np.tile(w1, (n_temp, 1))
        w2_temp_matrix = np.tile(w2_temp, (m, 1)).transpose()
        
        for vector1, vector2 in zip([x1, y1, z1], [x2_temp, y2_temp, z2_temp]):
            matrix1 = np.tile(vector1, (n_temp, 1))
            matrix2 = np.tile(vector2, (m, 1)).transpose()
            diff_mat += (matrix1 - matrix2)**2
            
        diff_mat = np.sqrt(diff_mat + (par_separation/10)**2)
        R += np.sum(1/diff_mat*w1_temp*w2_temp_matrix)

    return R

def Score_Overlap(x1, y1, z1, x2, y2, z2, w1, w2, par_separation):

    U12 = Potential_Overlap(x1, y1, z1, x2, y2, z2, w1, w2, par_separation)
    U11 = Potential_Overlap(x1, y1, z1, x1, y1, z1, w1, w1, par_separation)
    U22 = Potential_Overlap(x2, y2, z2, x2, y2, z2, w2, w2, par_separation)
    R = U12**2/U11/U22
    
    return R

def IsInBox(c, lowlim, highlim):
    c = np.array(c)
    lowlim = np.array(lowlim)
    highlim = np.array(highlim)
    diff_low = (c - lowlim) > 0
    diff_high = (c - highlim) < 0
    diff_sum = diff_low*diff_high
    if np.sum(diff_sum) < 3:
        return False
    else:
        return True
        
def SpuriousHalo_Mlim(mean_density, d, kmax):
    Mlim = 10.1*mean_density*d*kmax**(-2)
    return Mlim

def func_absdiff(y1, y2):
    return np.mean(np.abs(y1 - y2)**2)

def func_fitR(r, rho, rfactor):
    rho_r = rho*r**rfactor
    amp = np.mean(rho_r)
    return amp
    
def func_fitNFW_rho(r, rho):
    popt, pcov = curve_fit(func_logNFW, r, np.log10(rho))
    return popt, pcov

def func_fitPower_rho(r, rho):
    popt, pcov = curve_fit(func_logPowerLaw, r, np.log10(rho))
    return popt, pcov

def func_logPowerLaw(r, a, b):
    return np.log10(a) + b*np.log10(r)

def func_chisqr(obs, exp, error = 1, dof = 1):
    chisqr = np.sum( (obs - exp)**2 / error**2 )
    return chisqr/dof

def func_logPower32(log10r, a):
    return a + (-3/2)*log10r

def func_logPower127(log10r, a):
    return a + (-12/7)*log10r

def func_logNFW(log10r, rs, log10rho_s):
    return log10rho_s - log10r + np.log10(rs) - 2*np.log10(1 + 10**log10r/rs)

def func_logPowerPW(log10r, a, b):
    return a + b*log10r

def func_logBrokenNFW(log10r, log10rho_s, rs, a, b):
    return log10rho_s - a*log10r + a*np.log10(rs) - b*np.log10(1 + 10**log10r/rs)

def func_logBrokenPW(log10r, log10rho_s, log10rs, a, b):
    log10rho = np.zeros(len(log10r))
    small_r = log10r < log10rs
    log10rho[small_r] = log10rho_s - a*log10r[small_r] + a*log10rs
    log10rho[~small_r] = log10rho_s - b*log10r[~small_r] + b*log10rs
    
    return log10rho

def fit_all_functions(func_list, log_r, log_dens):
    """
    Fits density profiles by scanning the outer fitting radius (number of bins B)
    as described in Appendix A1 of the paper.
    It returns the parameters for the B that minimizes the reduced chi-squared.
    """
    
    # Initialize best-fit storage
    best_chi2 = np.ones(len(func_list)) * np.inf
    best_params = [None] * len(func_list)
    best_cov = [None] * len(func_list)
    available = np.zeros(len(func_list))
    
    # Minimum number of bins to fit (e.g., 20 as stated in paper)
    min_bins = 20
    total_bins = len(log_r)
    
    if total_bins < min_bins:
        # Not enough data points to even start scanning
        return available, best_params, best_cov

    # Scan B from min_bins to total_bins
    # This corresponds to varying the outer radius r_B
    for B in range(min_bins, total_bins + 1):
        
        # Slice data for the current radius r_B
        current_log_r = log_r[:B]
        current_log_dens = log_dens[:B]
        ini_dens = np.max(current_log_dens)
        
        for i, func in enumerate(func_list):
            try:
                p0 = None
                if (i == 0) or (i == 1):
                    p0 = [ini_dens]
                elif i == 2: # NFW
                    p0 = [np.max(10**current_log_r), ini_dens]
                elif i == 3: # Power Law Piecewise
                    p0 = [ini_dens, 0]
                elif i == 4: 
                    nfw_rs = np.max(10**current_log_r)/2 # Rough guess
                    nfw_rhos = ini_dens
                    p0 = [nfw_rs, nfw_rhos, 1, 2]
                elif i == 5: # Broken Power Law
                    p0 = [ini_dens, np.mean(current_log_r), -1.5, -5]

                param, cov = curve_fit(func, current_log_r, current_log_dens, p0=p0, maxfev=10000)
                
                fitted_y = func(current_log_r, *param)
                dof = len(current_log_r) - len(param)
                
                if dof > 0:
                    chi2 = np.sum((current_log_dens - fitted_y)**2)
                    red_chi2 = chi2 / dof
                    
                    if red_chi2 < best_chi2[i]:
                        best_chi2[i] = red_chi2
                        best_params[i] = param
                        best_cov[i] = cov
                        available[i] = 1
                        
            except RuntimeError:
                continue

    return available, best_params, best_cov

def get_powerlaw_envelope(r_data, q_data, d0):
            # Transform to linear space for fitting: Y = d0 - q
            Y = d0 - q_data
            # Filter for valid log inputs
            valid = (Y > 0) & (r_data > 0)
            if np.sum(valid) > 2:
                log_r = np.log(r_data[valid])
                log_Y = np.log(Y[valid])
                # Fit: log_Y = gamma * log_r + log_A
                coeffs = np.polyfit(log_r, log_Y, 1) 
                gamma = coeffs[0]
                A = np.exp(coeffs[1])
                return A, gamma
            return None, None

def find_min_peak(par_x0, par_y0, par_z0, par_linrho0, par_id_inhalo):

    par_x0_sorted = np.sort(np.unique(par_x0))
    spacing = par_x0_sorted[1] - par_x0_sorted[0]

    par_x0_size = int((np.max(par_x0) - np.min(par_x0))/spacing)
    par_y0_size = int((np.max(par_y0) - np.min(par_y0))/spacing)
    par_z0_size = int((np.max(par_z0) - np.min(par_z0))/spacing)
    n_grid = np.max([par_x0_size, par_y0_size, par_z0_size])+1
    grid = np.zeros((n_grid, n_grid, n_grid))
    id_grid = np.zeros((n_grid, n_grid, n_grid))
    
    x0_ongrid = ((par_x0 - np.min(par_x0))/spacing).astype(int)
    y0_ongrid = ((par_y0 - np.min(par_y0))/spacing).astype(int)
    z0_ongrid = ((par_z0 - np.min(par_z0))/spacing).astype(int)

    grid_coor = np.array(range(n_grid))
    x_grid_coord = grid_coor*spacing + np.min(par_x0)
    y_grid_coord = grid_coor*spacing + np.min(par_y0)
    z_grid_coord = grid_coor*spacing + np.min(par_z0)

    grid[x0_ongrid, y0_ongrid, z0_ongrid] = par_linrho0
    id_grid[x0_ongrid, y0_ongrid, z0_ongrid] = par_id_inhalo

    filter_size = np.array(range(2, int(n_grid/2)))
    npeak_filter = np.zeros(len(filter_size))   

    for filter_index, size_filter in enumerate(filter_size):
        coordinates = peak_local_max(grid, min_distance=size_filter)

        npeak_filter[filter_index] = coordinates.shape[0]

    npeak_gradient = np.gradient(npeak_filter)
    is_npeak_flat = np.sum(npeak_gradient == 0)

    if is_npeak_flat == 0:
        filter_size_final = np.max(filter_size)
    else:
        filter_size_final = np.min(filter_size[npeak_gradient == 0])

    peaks_final = peak_local_max(grid, min_distance=filter_size_final)

    return spacing, filter_size_final, filter_size, npeak_filter, grid, id_grid, peaks_final, \
            x_grid_coord, y_grid_coord, z_grid_coord

        
def set_grid(par_x0, par_y0, par_z0, par_linrho0, par_id_inhalo):
    par_x0_sorted = np.sort(np.unique(par_x0))
    spacing = np.max(par_x0_sorted[1:] - par_x0_sorted[:-1])
    par_x0_size = int(np.ceil((np.max(par_x0) - np.min(par_x0))/spacing + 1))
    par_y0_size = int(np.ceil((np.max(par_y0) - np.min(par_y0))/spacing + 1))
    par_z0_size = int(np.ceil((np.max(par_z0) - np.min(par_z0))/spacing + 1))
    # n_grid = np.max([par_x0_size, par_y0_size, par_z0_size])

    grid = np.zeros((par_x0_size, par_y0_size, par_z0_size))
    id_grid = np.zeros((par_x0_size, par_y0_size, par_z0_size), dtype = np.int64)
    sizes = [par_x0_size, par_y0_size, par_z0_size]

    x0_ongrid = ((par_x0 - np.min(par_x0))/spacing).astype(int)
    y0_ongrid = ((par_y0 - np.min(par_y0))/spacing).astype(int)
    z0_ongrid = ((par_z0 - np.min(par_z0))/spacing).astype(int)

    # grid_coor = np.array(range(n_grid))
    x_grid_coord = np.array(range(par_x0_size))*spacing + np.min(par_x0)
    y_grid_coord = np.array(range(par_y0_size))*spacing + np.min(par_y0)
    z_grid_coord = np.array(range(par_z0_size))*spacing + np.min(par_z0)

    grid[x0_ongrid, y0_ongrid, z0_ongrid] = par_linrho0
    id_grid[x0_ongrid, y0_ongrid, z0_ongrid] = par_id_inhalo

    return grid, id_grid, x_grid_coord, y_grid_coord, z_grid_coord, sizes

def next_step(x, y, z, grid, r):
    local_grid = grid[x-r:x+(r+1), y-r:y+(r+1), z-r:z+(r+1)]
    
    size_local = local_grid.shape[0]*local_grid.shape[1]*local_grid.shape[2]
    # print(local_grid.shape)
    
    if size_local < (2*r+1)**3:
        return 0.5, 0.5, 0.5
    # ix, iy, iz = np.argmax(local_grid, keepdims=True)
    ix, iy, iz = np.where(local_grid == np.max(local_grid))
    
    if len(ix) > 1:
        return ix[0], iy[0], iz[0]
    
    return ix, iy, iz

def find_max(start, grid, N_max, r):

    temp = start
    tem_list = np.zeros((3, N_max))
    for istep in range(N_max):
        # print(temp.shape)
        ix, iy, iz = next_step(temp[0], temp[1], temp[2], grid, r)
        
        if ix == 0.5:
            return np.array([np.nan, np.nan, np.nan]), tem_list
        
        # print(ix, iy, iz)
        temp = np.array([temp[0] + np.min(ix) - r, temp[1] + np.min(iy) - r, temp[2] + np.min(iz) - r])
        tem_list[0, istep] = temp[0]
        tem_list[1, istep] = temp[1]
        tem_list[2, istep] = temp[2]
        
        if (temp[0] > r)*(temp[1] > r)*(temp[2] > r) == 0:
            return np.array([np.nan, np.nan, np.nan]), tem_list
        elif (np.abs(grid.shape[0] - temp[0]) > r)*(np.abs(grid.shape[1] - temp[1]) > r)*(np.abs(grid.shape[2] - temp[2]) > r) == 0:
            return np.array([np.nan, np.nan, np.nan]), tem_list
        
        elif (np.min(np.abs(ix - r)) + np.min(np.abs(iy - r)) + np.min(np.abs(iz - r))) == 0:
            return temp, tem_list
            
    # return np.array([np.nan, np.nan, np.nan]), tem_list
    return temp, tem_list

def find_max_grid(grid, N_points, r, collapse_lim, N_max):
    final_list = np.zeros((3, N_points))
    for ip in range(N_points):

        rand_x = np.random.randint(low = r, high = grid.shape[0]-r)
        rand_y = np.random.randint(low = r, high = grid.shape[1]-r)
        rand_z = np.random.randint(low = r, high = grid.shape[2]-r)

        if grid[rand_x, rand_y, rand_z] > collapse_lim:
            point, traj = find_max(np.array([rand_x, rand_y, rand_z]), grid, N_max, r)

            final_list[0, ip] = point[0]
            final_list[1, ip] = point[1]
            final_list[2, ip] = point[2]
        else:
            final_list[:, ip] = np.nan
    return final_list


def func_show_long(long_dict):
    example_part = {key: long_dict[key] for key in list(long_dict.keys())[:10]}
    for key, value in example_part.items():
        print(f"Key: {key}, Value: {value}")

    return 0