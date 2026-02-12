import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import warnings; warnings.simplefilter('ignore')
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
import pandas as pd
from scipy.stats.stats import spearmanr
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

k_factor = (ks/kfs)**2
D_factor = (1-2/3*k_factor)*np.exp(-k_factor)
delta_sf = delta*D_factor
logdelta = np.log10(delta_sf)
ks_l = ks[logdelta > -2]
logdelta_l = logdelta[logdelta > -2]
logdelta2 = 2*logdelta_l
max_delta = np.argmax(logdelta2)
k_max = ks_l[max_delta]
print('k_max: ', k_max)

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
    
    
    File = Read_file_ini(path, 256, linearDensity)
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
    print(label, path)
    Par_Coor[label] = File[0]/Hconst
    Par_Mass[label] = File[1]*10**10/Hconst
    Par_ID[label] = File[2]
    # Par_Vel[label] = File[3]/Hconst

    if linearDensity is True:
        Par_LinearDensity[label] = File[4]
        Redshift[label] = File[5]
    else:
        Redshift[label] = File[4]




Par_FirstIn_T = {}
Par_SubIndex = {}
Indexer_Halo = {}
Halo_Dist = {}
Dist_Axis = {}
Halo_Density = {}

Par_FirstIn_T = {}
Par_FirstSubIndex = {}

SubPeak_Index = {}
SubPeak_Score = {}
SubHalo_ParCoor0 = {}
SubHalo_ParCoor = {}
SubHalo_LagCA = {}
LinRhoAllHaloPar = {}
Halo_Index = {}

Peaks_Properties = {}
CoorAllHaloPar0 = {}
Peak_IDs = {}
Semi_Selection = {}

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_min)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Peak_Data/AllPeak_{}.npy'.format(label)
try:
    Peak_data = np.load(path, allow_pickle = True)
    
    Peaks_Properties.update(Peak_data[0])
    CoorAllHaloPar0.update(Peak_data[1])
    # LinRhoAllHaloPar.update(Peak_data[2])
    # Peak_IDs.update(Peak_data[3])
    Semi_Selection.update(Peak_data[4])
    print(path)
except:
    print('No peak data')
    exit()

path = '/cosma7/data/dp004/dc-wang4/VVV_data/Particle_Data/AllParHistory_{}.npy'.format(label)
try:
    AllPar_data = np.load(path, allow_pickle = True)
    
    Par_FirstIn_T.update(AllPar_data[0])
    Par_FirstSubIndex.update(AllPar_data[1])
    print(path)
except:
    print('No ParHistory data')
    exit()


for T_min_temp in range(0, T_min+1):
    # for T_max_temp in range(200, T_min_temp-1, -1):
    label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min_temp, T_min_temp)

    path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/AllSubRef_{}.npy'.format(label)
    try:
        AllHalo_data = np.load(path, allow_pickle = True)
        
        # Indexer_Halo.update(AllHalo_data[0])
        Halo_Dist.update(AllHalo_data[1])
        Dist_Axis.update(AllHalo_data[2])
        Halo_Density.update(AllHalo_data[3])
        Par_SubIndex.update(AllHalo_data[4])
        Halo_Index.update(AllHalo_data[11])
        print(path)
    except:
        pass

    path = '/cosma7/data/dp004/dc-wang4/VVV_data/SubPeak_Data/AllSubPeak_{}.npy'.format(label)
    try:
        AllSub_data = np.load(path, allow_pickle = True)
        
        SubPeak_Index.update(AllSub_data[0])
        SubPeak_Score.update(AllSub_data[1])
        SubHalo_ParCoor0.update(AllSub_data[2])
        SubHalo_ParCoor.update(AllSub_data[3])
        SubHalo_LagCA.update(AllSub_data[4])
        # LinRhoAllHaloPar.update(AllSub_data[5])
        print(path)
    except:
        pass

        # path = '/cosma7/data/dp004/dc-wang4/VVV_data/Peak_Data/AllPeak_{}.npy'.format(label)
        # try:
        #     Peak_data = np.load(path, allow_pickle = True)
            
        #     Peaks_Properties.update(Peak_data[0])
        #     CoorAllHaloPar0.update(Peak_data[1])
        #     Peak_IDs.update(Peak_data[3])
        #     Semi_Selection.update(Peak_data[4])
        #     print(path)
        # except:
        #     continue
           
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
        # SubHalo_Mass[label] = File[5]*10**10/Hconst
        Halo_FirSub[label] = File[6]
        Halo_NSub[label] = File[7]
        # SubHalo_GrN[label] = File[8]
        Halo_Offset[label] = File[9]
        Halo_ParN[label] = File[10]
        SubHalo_Offset[label] = File[11]
        SubHalo_ParN[label] = File[12]
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
        T_max = T - 1
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

box_size = Box_size[label]
Chosen_Peak = {}
Halo_CausticsR = {}
LocalRho_Halo = {}
LocalRho_Semi = {}

RhoFitParam_Halo = {}
RhoFitParam_Semi = {}
RhoFitChi2_Halo = {}
RhoFitChi2_Semi = {}
Spear_Coef = {}
HaloMass_Data = {}
LocalRho_Various = {}
Peak_Density = {}
Peak_Dist_toCentre = {}
Ell_Comp = {}
Ell_LargestR = {}
Selections = {}


npeakmax = 20
Ndistbin = 50
r_high_lim = halor
N_nfw = 51
N_rbin = N_nfw
Na = 10
NbinH = 70
dens_factorH = 100
Nrand = int(1e3)
filter_window = 7
filter_order = 3
rfactor = 12/7

label0 = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, -1)
a_ini = 1/(Redshift[label0]+1)
# Par_x0 = Par_Coor[label0][:, 0]
# Par_y0 = Par_Coor[label0][:, 1]
# Par_z0 = Par_Coor[label0][:, 2]
# Par_m0 = Par_Mass[label0]
# Par_id0 = Par_ID[label0]
# Par_linrho = Par_LinearDensity[label0]
index_parent0 = pd.Index(Par_ID[label0])

high_selection = Par_Mass[label0] == Par_Mass[label0][0]
Par_x0_m = Par_Coor[label0][:, 0][high_selection]
# Par_y0_m = Par_Coor[label0][:, 1][high_selection]
# Par_z0_m = Par_Coor[label0][:, 2][high_selection]
range_x = np.max(Par_x0_m) - np.min(Par_x0_m)
# range_y = np.max(Par_y0_m) - np.min(Par_y0_m)
# range_z = np.max(Par_z0_m) - np.min(Par_z0_m)

# par_separation = np.mean([range_x, range_y, range_z])/(Nh*20)
par_separation = range_x/(Nh*20)

mean_density0 = np.sum(Par_Mass[label0])/542.16**3
halo_mlim = SpuriousHalo_Mlim(mean_density0, par_separation, k_max*Hconst)

print(halo_mlim)

nrows = 3
ncols = 6
irow_dens_halo = (0, 0)
irow_firstTdens = (0, 1)
irow_surr = (0, 2)
irow_spu = (0, 3)
irow_nonsemi = (0, 4)
irow_noavailable = (0, 5)

irow_ini = (1, 0)
irow_lin = (1, 1)
irow_lin_peak = (1, 2)
irow_lin_semi = (1, 3)
irow_dens_peak = (1, 4)
irow_subpeak = (1, 5)

irow_dens = (2, 0)
irow_caustics = (2, 1)
irow_dens_centre = (2, 2)
irow_dens_ava_centre = (2, 3)
irow_dens_fit = (2, 5)

a_list = np.zeros(T_max - T_min + 1)
T_index = 0
func_list = [func_logPower32, func_logPower127, func_logNFW, func_logPowerPW, func_logBrokenNFW, func_logBrokenPW]
func_color_list = ['black', 'lightskyblue', 'mediumpurple', 'mediumseagreen', 'orange', 'teal']

for T in range(T_min, T_max+1):
# for T in [60]:
    plt.clf()
    fig, axs = plt.subplots(figsize = (ncols*2+4, nrows*2+4),  #sharex = 'row', sharey = 'row',
                            nrows = nrows, ncols = ncols, 
                            dpi = 200)
    # axs = axs.flat
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
    print(label)
    a_list[T_index] = 1/(Redshift[label]+1)
    a_temp = a_list[T_index]
    print(a_temp)

    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)

    if len(Halo_Mass_200[label]) == 0:
        print('No Halo at ', label)
        continue
    try:
        halo_index = int(Halo_Index[label])
    except:
        print('No progenitor at ', label)
        continue

    if len(Halo_NSub[label]) == 0:
        continue
    halo_x = Halo_Coor[label][halo_index, 0]*a_temp
    halo_y = Halo_Coor[label][halo_index, 1]*a_temp
    halo_z = Halo_Coor[label][halo_index, 2]*a_temp
    halo_r = Halo_R_200[label][halo_index]*a_temp
    halo_pos = np.array([halo_x, halo_y, halo_z])
    halo_firstsub = Halo_FirSub[label][halo_index]
    halo_nsub = Halo_NSub[label][halo_index]
    halo_off = Halo_Offset[label][halo_index]
    halo_np = Halo_ParN[label][halo_index]
    HaloMass_Data[label] = Halo_Mass_200[label][halo_index]
    if halo_nsub == 0:
        print ('no subhalo in:', label)
        continue
        
    # Par_x_cov = Par_Coor[label][:, 0]
    # Par_x = Par_x_cov*a_temp
    # Par_y_cov = Par_Coor[label][:, 1]
    # Par_y = Par_y_cov*a_temp
    # Par_z_cov = Par_Coor[label][:, 2]
    # Par_z = Par_z_cov*a_temp
    
    # Par_vx = Par_Vel[label][:, 0]*a_temp
    # Par_vy = Par_Vel[label][:, 1]*a_temp
    # Par_vz = Par_Vel[label][:, 2]*a_temp
    # Par_id = Par_ID[label]
    # Par_m = Par_Mass[label]

    par_x_inhalo = Par_Coor[label][int(halo_off):int(halo_off+halo_np), 0]*a_temp
    par_y_inhalo = Par_Coor[label][int(halo_off):int(halo_off+halo_np), 1]*a_temp
    par_z_inhalo = Par_Coor[label][int(halo_off):int(halo_off+halo_np), 2]*a_temp
    par_m_inhalo = Par_Mass[label][int(halo_off):int(halo_off+halo_np)]
    par_id_inhalo = Par_ID[label][int(halo_off):int(halo_off+halo_np)]
    par_vx_inhalo = Par_Vel[label][int(halo_off):int(halo_off+halo_np), 0]*a_temp
    par_vy_inhalo = Par_Vel[label][int(halo_off):int(halo_off+halo_np), 1]*a_temp
    par_vz_inhalo = Par_Vel[label][int(halo_off):int(halo_off+halo_np), 2]*a_temp
    
    mean_density = np.sum(Par_Mass[label])/(Box_size[label]*a_temp)**3
    mass_crit_ratio = (mean_density/Par_Mass[label][0])
    
    ipeak_index = []
    ipeak_dist = []
    npeakmax = np.min([Peaks_Properties[label], npeakmax])
    print('Number of peaks: ', Peaks_Properties[label])
    normal_peak = colors.Normalize(vmin = 0, vmax = npeakmax)
    for ipeak in range(npeakmax):
        label_peak = '{}P{}'.format(label, ipeak)
        try:
            [centre_maxlinrho_temp, peak_maxlinrho_temp, 
            dist_near_peak, linrho_near_peak, 
            quad_coef, log_coefs, ell_fit_coef, 
            peak_r, ellipticity, prolateness, 
            eigenvals, eigenvecs, ell_fit_matrix,
            id_temp, x_temp, y_temp, z_temp] = Peaks_Properties[label_peak]


            [semi_selection, par_x_semi, par_y_semi, par_z_semi] = Semi_Selection[label_peak]

        except:
            print('Peak semi selection not exist!')
            continue
        try:
            dist_halo = Halo_Dist[label]
            [dist_bin, d_dist, dist_bin_mid] = Dist_Axis[label]
        except:
            print('Halo dist data not exist!')
            continue
            
        print(label_peak, np.sum(semi_selection)/len(semi_selection))

        halo_count, subbin_edges = np.histogram(dist_halo[semi_selection], bins = dist_bin)
        density_peak = halo_count/d_dist/dist_bin_mid**2/4/np.pi*Par_Mass[label][0]

        dist_peak_to_halo = np.sqrt(np.sum((np.array(peak_maxlinrho_temp)*a_temp - halo_pos)**2))
        Peak_Density[label_peak] = density_peak
        Peak_Dist_toCentre[label_peak] = dist_peak_to_halo
        x_axis = dist_bin_mid
        y_axis = (density_peak)*(x_axis)**(rfactor)
        axs[irow_dens_peak].plot(x_axis, y_axis, color = cm.rainbow(normal_peak(ipeak)), linewidth = 2, 
                                            alpha = 1, linestyle = 'solid')

        ipeak_dist.append(np.mean(np.log10(y_axis[:5])))
        ipeak_index.append(ipeak)
        

    if len(ipeak_dist) == 0:
        print('no peak available', label)
        exit()
        
    ipeak_dist = np.array(ipeak_dist)
    ipeak_index = np.array(ipeak_index)
    chosen_peak_index = ipeak_index[int(np.argmax(ipeak_dist))]
    
    
    label_peak = '{}P{}'.format(label, chosen_peak_index)
    print('choosing {}th peak'.format(chosen_peak_index-1), label_peak)
    
    # [centre_maxlinrho, peak_maxlinrho, 
    # dist_near_peak, linrho_near_peak, 
    # quad_coef, log_coefs, ell_fit_coef, 
    # peak_r, ellipticity, prolateness, 
    # id_near_maxlinrho, x_near_maxlinrho, y_near_maxlinrho, z_near_maxlinrho] 
    
    [centre_maxlinrho, peak_maxlinrho, 
    dist_near_peak, linrho_near_peak, 
    quad_coef, log_coefs, ell_fit_coef, 
    peak_r, ellipticity, prolateness, 
    eigenvals, eigenvecs, ell_fit_matrix,
    id_near_maxlinrho, x_near_maxlinrho, y_near_maxlinrho, z_near_maxlinrho] = Peaks_Properties[label_peak]

    func_scatter(axs[irow_lin_peak], dist_near_peak, linrho_near_peak, Nrand, 1, 0.5, c = 'grey')
    
    Chosen_Peak[label] = [chosen_peak_index, Peaks_Properties[label_peak]]

    par_x_nearmax_c = x_near_maxlinrho - centre_maxlinrho[0]
    par_y_nearmax_c = y_near_maxlinrho - centre_maxlinrho[1]
    par_z_nearmax_c = z_near_maxlinrho - centre_maxlinrho[2]
    dist_near_peak_ell = np.sqrt(func_FitQua_ellipsoid([par_x_nearmax_c, par_y_nearmax_c, par_z_nearmax_c, 0], 
                                        *ell_fit_coef)/quad_coef[0])
    func_scatter(axs[irow_lin_peak], dist_near_peak_ell, linrho_near_peak, Nrand, 1, 0.5, c = 'blue')
    
    index_bin_min, dist_bin_min, linrho_bin_min, \
    index_bin_max, dist_bin_max, linrho_bin_max = func_minmax(dist_near_peak_ell, linrho_near_peak, Ndistbin)
    
    axs[irow_lin_peak].plot(dist_bin_max, linrho_bin_max, color = 'blue', linewidth = 2)
    dist_axis = np.linspace(np.min(dist_near_peak), np.max(dist_near_peak), Ndistbin)
    quad_coef_upper = [0, 0]
    quad_coef_upper[0], quad_coef_upper[1] = func_Fitqualin(dist_bin_max, linrho_bin_max)
    quad_fit_upper = np.polyval(quad_coef_upper, dist_axis**2)
    axs[irow_lin_peak].plot(dist_axis, quad_fit_upper, color = 'red', linewidth = 2)
    
    [par_x_inhalo0, par_y_inhalo0, par_z_inhalo0] = CoorAllHaloPar0[label]
    par_inhalo_selection = index_parent0.get_indexer(par_id_inhalo)
    par_inhalo_selection = par_inhalo_selection[par_inhalo_selection >= 0]
    par_linrho_inhalo0 = Par_LinearDensity[label0][par_inhalo_selection]

    if len(par_inhalo_selection) < 100:
        continue
    
    par_dist_inhalo0 = distance(par_x_inhalo0, par_y_inhalo0, par_z_inhalo0, centre_maxlinrho)
    
    par_x_inhalo0_c = par_x_inhalo0 - centre_maxlinrho[0]
    par_y_inhalo0_c = par_y_inhalo0 - centre_maxlinrho[1]
    par_z_inhalo0_c = par_z_inhalo0 - centre_maxlinrho[2]
    dist_inhalo0_ell = np.sqrt(func_FitQua_ellipsoid([par_x_inhalo0_c, par_y_inhalo0_c, par_z_inhalo0_c, 0], 
                                        *ell_fit_coef)/quad_coef[0])
    
    func_scatter(axs[irow_lin], par_dist_inhalo0, par_linrho_inhalo0, Nrand, 0.1, 0.5, c = 'grey')
    # func_scatter(axs[irow_lin_semi], par_dist_inhalo0, par_linrho_inhalo0, Nrand, 0.1, 0.5, c = 'grey')
    func_scatter(axs[irow_ini], par_x_inhalo0, par_y_inhalo0, Nrand, 0.1, 0.1, c = 'grey')
    func_scatter(axs[irow_lin], dist_inhalo0_ell, par_linrho_inhalo0, Nrand, 0.1, 0.1, c = 'blue')
    func_scatter(axs[irow_lin_semi], dist_inhalo0_ell, par_linrho_inhalo0, Nrand, 0.1, 0.1, c = 'blue')
    
    dist_axis_halo = np.linspace(np.min(dist_inhalo0_ell), np.max(dist_inhalo0_ell), Ndistbin)
    # quad_fit_upper_halo_axis = np.polyval(quad_coef_upper, dist_axis_halo**2)
    axs[irow_lin].plot(dist_axis_halo, quad_fit_upper_halo_axis, color = 'red', linewidth = 2)
    
    r_fit_max = 0.4 * peak_r
    fit_mask = (dist_near_peak_ell > 0) & (dist_near_peak_ell <= r_fit_max)
    
    dist_for_fit = dist_near_peak_ell[fit_mask]
    linrho_for_fit = linrho_near_peak[fit_mask]

    Nbins_method = 50
    bin_edges = np.linspace(np.min(dist_for_fit), np.max(dist_for_fit), Nbins_method + 1)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    bin_indices = np.digitize(dist_for_fit, bin_edges)

    q01_env = []  # Lower envelope 
    q99_env = []  # Upper envelope 
    r_env = []    

    for i in range(1, Nbins_method + 1):
        in_bin = bin_indices == i
        if np.sum(in_bin) > 10:  
            deltas_in_bin = linrho_for_fit[in_bin]
            q01_env.append(np.quantile(deltas_in_bin, 0.01))
            q99_env.append(np.quantile(deltas_in_bin, 0.99))
            r_env.append(bin_centers[i-1])
            
    r_env = np.array(r_env)
    q01_env = np.array(q01_env)
    q99_env = np.array(q99_env)
    
    delta0 = np.max(linrho_near_peak)

    if len(r_env) > 2:
        A_up, gamma_up = get_powerlaw_envelope(r_env, q99_env, delta0)
        A_lo, gamma_lo = get_powerlaw_envelope(r_env, q01_env, delta0)

        fit_upper = delta0 - A_up * dist_inhalo0_ell**gamma_up
        fit_lower = delta0 - A_lo * dist_inhalo0_ell**gamma_lo

        semi_selection = (par_linrho_inhalo0 < fit_upper) & (par_linrho_inhalo0 > fit_lower)
        
    else:
        print("Not enough points to fit envelopes")
        semi_selection = np.zeros(len(par_linrho_inhalo0), dtype=bool)

    # quad_fit_upper_halo = np.polyval(quad_coef_upper, dist_inhalo0_ell**2)
    # quad_coef_lower = [0, 0]
    # quad_coef_lower[0], quad_coef_lower[1] = func_Fitqualin(dist_bin_min, linrho_bin_min)
    # quad_fit_lower_halo = np.polyval(quad_coef_lower, dist_inhalo0_ell**2)

    # semi_selection = (par_linrho_inhalo0 < quad_fit_upper_halo) & (par_linrho_inhalo0 > quad_fit_lower_halo)
    # semi_selection = par_linrho_inhalo0 < quad_fit_upper_halo
    func_scatter(axs[irow_lin_semi], dist_inhalo0_ell, quad_fit_upper_halo, Nrand, 0.1, 0.5, c = 'red')
    func_scatter(axs[irow_lin_semi], dist_inhalo0_ell[semi_selection], 
                 par_linrho_inhalo0[semi_selection], Nrand, 0.1, 0.5, c = 'salmon')
    
    dens_map, xedges, yedges = np.histogram2d(par_x_inhalo, par_y_inhalo, bins=[NbinH, NbinH],
                                             range = [[halo_x - halo_r, halo_x + halo_r], 
                                                     [halo_y - halo_r, halo_y + halo_r]]) 
    
    # par_selection = func_SphSelect(Par_x, Par_y, Par_z, halo_r*halor,
    #                                 halo_x, halo_y, halo_z)
    # extent = [halo_x - halo_r, halo_x + halo_r, halo_y - halo_r, halo_y + halo_r]
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    axs[irow_dens_halo].imshow(np.transpose(np.log10(dens_map+dens_factorH)), 
                                extent=extent, origin='lower',
                                cmap=cm.twilight_shifted, interpolation="spline16")
    axs[irow_dens_halo].scatter(peak_maxlinrho[0]*a_temp, peak_maxlinrho[1]*a_temp, 
                                alpha = 0.8, 
                                marker = 'x', s = 100, c = 'red')
    axs[irow_dens_halo].scatter(halo_x, halo_y, 
                                alpha = 0.8, 
                                marker = '+', s = 100, c = 'black')
    
    dens_map_semi, xedges, yedges = np.histogram2d(par_x_inhalo[~semi_selection], 
                                                  par_y_inhalo[~semi_selection], 
                                              bins=[NbinH, NbinH],
                                              range = [[halo_x - halo_r, halo_x + halo_r], 
                                                       [halo_y - halo_r, halo_y + halo_r]]) 
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    axs[irow_nonsemi].imshow(np.transpose(np.log10(dens_map_semi+np.mean(dens_map_semi))), 
                                    extent=extent, origin='lower',
                                    cmap=cm.twilight_shifted, interpolation="spline16")
    
    
    [dist_bin_mid, density_halo] = Halo_Density[label]
    # [halo_thetas, halo_phis, dist_bin_mid_ag, halo_dens_r_ag, 
    #  halo_dens_r_err, halo_dens_r_agbin] = LocalRhoAg_Halo[label]
    print('Density:', len(dist_bin_mid), len(density_halo))
    
    dist_bin_mid = dist_bin_mid*a_temp
    # dist_bin_mid_ag = dist_bin_mid_ag*a_temp
    density_halo = density_halo/a_temp**3
    # halo_dens_r_ag = halo_dens_r_ag/a_temp**3
    
    x_axis = dist_bin_mid
    y_axis = (density_halo)*(x_axis)**(rfactor)
    axs[irow_dens].plot(x_axis, y_axis, color = 'black', 
                        linewidth = 2, alpha = 1, linestyle = 'solid')
    
    dlogrho_dlogr = func_HaloCaustics(dist_bin_mid, density_halo, 15, 4)

    is_local_min = (dlogrho_dlogr[1:-1] < dlogrho_dlogr[0:-2]) & \
                   (dlogrho_dlogr[1:-1] < dlogrho_dlogr[2:])
    
    local_min_indices = np.where(is_local_min)[0] + 1
    if len(local_min_indices) > 0:
        r_caustics = dist_bin_mid[local_min_indices[0]]
    else:
        r_caustics = dist_bin_mid[np.argmin(dlogrho_dlogr)]

    # r_caustics = dist_bin_mid[np.argmin(dlogrho_dlogr)]
    Halo_CausticsR[label] = r_caustics
    axs[irow_caustics].plot(dist_bin_mid, dlogrho_dlogr, color = 'black', 
                            linewidth = 2, alpha = 1, linestyle = 'solid')
    axs[irow_caustics].axvline(x = r_caustics, 
                            linewidth = 2, alpha = 1, linestyle = 'dashed')
    axs[irow_dens].axvline(x = r_caustics, 
                            linewidth = 2, alpha = 1, linestyle = 'dashed')
    axs[irow_dens].axvline(x = halo_r, color = 'red', 
                            linewidth = 2, alpha = 1, linestyle = 'dashed')

    Halo_CausticsR[label] = r_caustics
    
    
    log_r_halo = np.log10(dist_bin_mid[dist_bin_mid < r_caustics])
    log_dens_halo = np.log10(density_halo[dist_bin_mid < r_caustics])
    x_axis = 10**log_r_halo
    y_axis = 10**log_dens_halo*(x_axis)**(rfactor)
    axs[irow_dens_centre].plot(x_axis, y_axis, color = 'grey', linewidth = 2, alpha = 1, linestyle = 'solid')
    axs[irow_dens_ava_centre].plot(x_axis, y_axis, color = 'grey', linewidth = 2, alpha = 0.4, linestyle = 'solid')
    
    # axs[irow_dens_centre].axvline(x = halo_r, color = 'red', 
    #                             linewidth = 2, alpha = 1, linestyle = 'dashed')
    if len(log_r_halo) == 0:
        print('Pass density: ', label)
        continue
    LocalRho_Halo[label] = [log_r_halo, log_dens_halo]
    available_halo, param_list_halo, cov_list_halo = fit_all_functions(func_list, log_r_halo, log_dens_halo )
    
    RhoFitParam_Halo[label] = [available_halo, param_list_halo, cov_list_halo]

    fit_chi2_halo = np.zeros(len(func_list))
    index = 0
    for i, func in enumerate(func_list):
        if available_halo[i] == 1:
            param = param_list_halo[index]
            cov = cov_list_halo[index]
            fitted_y = func(log_r_halo, *param)
            x_axis = 10**log_r_halo
            y_axis = 10**fitted_y*(x_axis)**(rfactor)
            axs[irow_dens_centre].plot(x_axis, y_axis, color = func_color_list[i], 
                                linewidth = 2, alpha = 1, linestyle = 'solid')
            fit_chi2_halo[i] = func_chisqr(obs = log_dens_halo, exp = fitted_y, 
                                           dof = len(x_axis) - len(param))
            # chi2, fit_chi2_halo[i] = stats.chisquare(f_obs=log_dens_halo, f_exp=fitted_y, 
            #                                   ddof=len(x_axis) - len(param))
            index += 1
    
    axs[irow_dens_fit].plot(np.arange(len(func_list))[available_halo == 1], 
                            fit_chi2_halo[available_halo == 1], color = 'black', marker = 'x')
    RhoFitChi2_Halo[label] = fit_chi2_halo

    par_first_inhalo_T = Par_FirstIn_T[label]
    par_firstin_subindex = Par_FirstSubIndex[label]
    print(len(par_first_inhalo_T), len(par_firstin_subindex))
    
    T_map, xedges, yedges = np.histogram2d(par_x_inhalo, par_y_inhalo, bins=[NbinH, NbinH], 
                                               weights = par_first_inhalo_T - np.min(par_first_inhalo_T),
                                               range = [[halo_x - halo_r, halo_x + halo_r], 
                                               [halo_y - halo_r, halo_y + halo_r]]) 
    
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    densT_map = T_map/dens_map
    for index in np.argwhere(np.isnan(densT_map)):
        m = index[0]
        n = index[1]
        densT_map[m, n] = 0
    axs[irow_firstTdens].imshow(np.transpose(np.log10(densT_map + np.mean(densT_map))), 
                                extent=extent, origin='lower',
                                cmap=cm.twilight_shifted, interpolation="spline16")
    
    par_selection_surr = func_SphSelect(Par_Coor[label][:, 0], 
                                        Par_Coor[label][:, 1], 
                                        Par_Coor[label][:, 2], halo_r*halor, #peak_x, peak_y, peak_z)
                                    halo_x, halo_y, halo_z)
    
    par_x_inhalo_surr = Par_Coor[label][:, 0][par_selection_surr]
    par_y_inhalo_surr = Par_Coor[label][:, 1][par_selection_surr]
    par_z_inhalo_surr = Par_Coor[label][:, 2][par_selection_surr]

    par_id_inhalo_surr = Par_ID[label][par_selection_surr]
    par_m_inhalo_surr = Par_Mass[label][par_selection_surr]
    dens_map, xedges, yedges = np.histogram2d(par_x_inhalo_surr, par_y_inhalo_surr, bins=[NbinH, NbinH]) 
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    axs[irow_surr].imshow(np.transpose(np.log10(dens_map+dens_factorH)), 
                                extent=extent, origin='lower',
                                cmap=cm.twilight_shifted, interpolation="spline16")
    
    
    par_isspurious = np.ones(len(par_first_inhalo_T))
    
    for T_min_temp in range(int(np.min(par_first_inhalo_T)), int(np.max(par_first_inhalo_T)+1)):
        par_firstin_subindex_temp = np.zeros(len(par_first_inhalo_T)) - 1
        
        label_temp = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_min_temp)
        # halo_nsub = Halo_NSub[label][halo_index]
        [halo_nsub, par_sub_index] = Par_SubIndex[label_temp]
        Ntry_sub = np.min([halo_nsub, 20])
        sub_with_peak_list = []
        spu_list = []
        
        par_firstT_selection = par_first_inhalo_T == T_min_temp
        if np.sum(par_firstT_selection) == 0:
            continue
        par_firstin_subindex_temp[par_firstT_selection] = par_firstin_subindex[par_firstT_selection]
        
        for isub in range(Ntry_sub):
            label_sub = 'L{}HaloID{}k{}Nh{}S{}SH{}'.format(L, HaloID, kfs, Nh, T_min_temp, isub)
            peak_overlap = SubPeak_Score[label_sub]
            subpeak_index = SubPeak_Index[label_sub]
            npeak_insub = len(subpeak_index)
            if isub == 0:
                mean_score = np.median(peak_overlap)
            
            if (npeak_insub > 0) and (np.max(peak_overlap) > mean_score):
            # if (npeak_insub > 0):
                sub_with_peak_list.append(isub)
            elif (npeak_insub > 0) and (np.max(peak_overlap) <= mean_score):
                spu_list.append(isub)
        index_subpeak_parent = pd.Index(sub_with_peak_list)
        par_subpeak_indexer = index_subpeak_parent.get_indexer(par_firstin_subindex_temp)
        par_subpeak_selection = par_subpeak_indexer >= 0

        index_spu_parent = pd.Index(spu_list)
        par_spu_indexer = index_spu_parent.get_indexer(par_firstin_subindex_temp)
        par_spu_selection = par_spu_indexer >= 0
        
        par_isspurious[par_firstT_selection*par_spu_selection] = 2
        par_isspurious[par_firstT_selection*par_subpeak_selection] = 0
        
    print('Fraction confirmed spurious:', np.sum(par_isspurious == 2)/len(par_isspurious))
    print('Fraction confirmed major subhalo:', np.sum(par_isspurious == 0)/len(par_isspurious))
    print('Fraction smooth:', np.sum(par_isspurious == 1)/len(par_isspurious))
        
    
    halo_spurious_selection = par_isspurious > 1
    dens_map_spu, xedges, yedges = np.histogram2d(par_x_inhalo[halo_spurious_selection], 
                                                  par_y_inhalo[halo_spurious_selection], 
                                              bins=[NbinH, NbinH],
                                              range = [[halo_x - halo_r, halo_x + halo_r], 
                                                       [halo_y - halo_r, halo_y + halo_r]]) 
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    axs[irow_spu].imshow(np.transpose(np.log10(dens_map_spu+np.mean(dens_map_spu))), 
                                    extent=extent, origin='lower',
                                    cmap=cm.twilight_shifted, interpolation="spline16")
    

    
    
    [halo_nsub, par_sub_index] = Par_SubIndex[label]
    normal = colors.Normalize(vmin = 0, vmax = Ntry_sub)
    for isub in range(Ntry_sub):
        label_sub = 'L{}HaloID{}k{}Nh{}S{}SH{}'.format(L, HaloID, kfs, Nh, T, isub)
        sub_index = halo_firstsub + isub
        par_insub_off = SubHalo_Offset[label][int(sub_index)]
        par_insub_np = SubHalo_ParN[label][int(sub_index)]
        par_off_low = int(par_insub_off - halo_off)
        par_off_high = int(par_insub_off - halo_off + par_insub_np)
        par_linrho_insub_temp = par_linrho_inhalo0[par_off_low:par_off_high]
            
        subpeak_index = SubPeak_Index[label_sub]
        peak_overlap = SubPeak_Score[label_sub]
        npeak_insub = len(subpeak_index)
        
        if (npeak_insub == 0):
            color_sub = 'black'
        
        else:
            
            peak_insub_index = np.argmax(peak_overlap)
            ipeak_insub = subpeak_index[peak_insub_index]
            score_insub = peak_overlap[peak_insub_index]
            
            if score_insub < mean_score:
                color_sub = 'black'
            else:
                color_sub = 'red'
            
            axs[irow_subpeak].scatter(ipeak_insub, score_insub,
                    color = cm.rainbow(normal(isub)), marker = 'x', s = 100)
            
            axs[irow_subpeak].scatter(subpeak_index, peak_overlap,
                    color = cm.rainbow(normal(isub)), marker = '+', s = 50)
        
        
        sub_ca = SubHalo_LagCA[label_sub]
        [par_x_insub0_temp, par_y_insub0_temp, par_z_insub0_temp] = SubHalo_ParCoor0[label_sub]
        [par_x_insub_temp, par_y_insub_temp, par_z_insub_temp] = SubHalo_ParCoor[label_sub]
        
        func_scatter(axs[irow_ini], par_x_insub0_temp, par_y_insub0_temp, Nrand, 0.1, 0.1, c = color_sub)
        dist0_sub_temp = distance(par_x_insub0_temp, par_y_insub0_temp, par_z_insub0_temp, centre_maxlinrho)
        func_scatter(axs[irow_lin], dist0_sub_temp, par_linrho_insub_temp, Nrand, 0.1, 0.5, c = color_sub)
        
    
    
    ava_selection = semi_selection * (~halo_spurious_selection)
    print('semi fraction:', np.sum(semi_selection)/len(semi_selection))
    print('available fraction:', np.sum(ava_selection)/len(ava_selection))
    if np.sum(ava_selection)/len(ava_selection) == 0:
        continue
    total_selection = ~ava_selection
    dens_map_available, xedges, yedges = np.histogram2d(par_x_inhalo[ava_selection], 
                                                  par_y_inhalo[ava_selection], 
                                              bins=[NbinH, NbinH],
                                              range = [[halo_x - halo_r, halo_x + halo_r], 
                                                       [halo_y - halo_r, halo_y + halo_r]]) 
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    axs[irow_noavailable].imshow(np.transpose(np.log10(dens_map_available+np.mean(dens_map_available))), 
                                    extent=extent, origin='lower',
                                    cmap=cm.twilight_shifted, interpolation="spline16")
    
    dist_halo = Halo_Dist[label]
    dist_halo = dist_halo*a_temp
    [dist_bin, d_dist, dist_bin_mid_temp] = Dist_Axis[label]
    dist_bin = dist_bin*a_temp
    d_dist = d_dist*a_temp
    color_sel_list = ['grey', 'brown', 'red', 'green', 'teal', 'orange', 'salmon']
    
    Selections[label] = [(~semi_selection)+semi_selection, ~semi_selection, semi_selection, halo_spurious_selection, total_selection, ava_selection]
    dens_all = []
    for i, selection in enumerate([(~semi_selection)+semi_selection, ~semi_selection, semi_selection, halo_spurious_selection, total_selection, ava_selection]):
        
        halo_count, subbin_edges = np.histogram(dist_halo[selection], bins = dist_bin)
        density_temp = halo_count/d_dist/dist_bin_mid**2/4/np.pi*par_m_inhalo[0]
        # Halo_Density[label] = [dist_bin_mid_temp, density_halo]
        x_axis = dist_bin_mid
        y_axis = (density_temp)*(x_axis)**(rfactor)
        axs[irow_dens].plot(x_axis, y_axis, color = color_sel_list[i], linewidth = 2, 
                            alpha = 0.5, linestyle = 'dashed')
        dens_all.append(density_temp)
    LocalRho_Various[label] = dens_all
        
    dist_ava = dist_halo[ava_selection]
    spear = spearmanr(dist_inhalo0_ell, dist_halo)[0]
    spear_ava = spearmanr(dist_inhalo0_ell[ava_selection], dist_ava)[0]
    halo_count_ava, subbin_edges = np.histogram(dist_ava, bins = dist_bin)
    density_ava = halo_count_ava/d_dist/dist_bin_mid**2/4/np.pi*par_m_inhalo[0]
    dist_bin_mid_ava = dist_bin_mid[density_ava > 0]
    density_ava = density_ava[density_ava > 0]
    Spear_Coef[label] = [spear, spear_ava]
        
    log_r_ava = np.log10(dist_bin_mid_ava[dist_bin_mid_ava < r_caustics])
    log_dens_ava = np.log10(density_ava[dist_bin_mid_ava < r_caustics])
    if np.sum(log_dens_ava > 0) > 0:
        x_axis = 10**log_r_ava
        y_axis = 10**log_dens_ava*(x_axis)**(rfactor)
        axs[irow_dens_ava_centre].plot(x_axis, y_axis, color = 'grey', 
                                    linewidth = 2, alpha = 1, linestyle = 'solid')
        LocalRho_Semi[label] = [log_r_ava, log_dens_ava]

        available_ava, param_list_ava, cov_list_ava = fit_all_functions(func_list, log_r_ava, log_dens_ava )
        fit_chi2_ava = np.zeros(len(func_list))
        print(len(param_list_ava))
    
        RhoFitParam_Semi[label] = [available_ava, param_list_ava, cov_list_ava]

        index = 0
        for i, func in enumerate(func_list):
            if available_ava[i] == 1:
                param = param_list_ava[index]
                cov = cov_list_ava[index]
                fitted_y = func(log_r_ava, *param)
                x_axis = 10**log_r_ava
                y_axis = 10**fitted_y*(x_axis)**(rfactor)
                axs[irow_dens_ava_centre].plot(x_axis, y_axis, color = func_color_list[i], 
                                    linewidth = 2, alpha = 1, linestyle = 'dashed')
                fit_chi2_ava[i] = func_chisqr(obs = log_dens_ava, exp = fitted_y, 
                                            dof = len(x_axis) - len(param))
                # chi2, fit_chi2_ava[i] = stats.chisquare(f_obs=log_dens_ava, f_exp=fitted_y, 
                #                                   ddof=len(x_axis) - len(param))
                index += 1
        
        axs[irow_dens_fit].plot(np.arange(len(func_list)), fit_chi2_ava, color = 'red', marker = '+')
        RhoFitChi2_Semi[label] = fit_chi2_ava
    else:
        print('No available particles!')

    for irow in [irow_dens, irow_dens_centre, irow_dens_ava_centre, irow_dens_peak]:
        axs[irow].set_xscale("log")
        axs[irow].set_yscale("log")
        axs[irow].set_xlabel(r'$r[Mpc]$')
    for irow in [irow_caustics]:
        axs[irow].set_xscale("log")
        axs[irow].set_xlabel(r'$r[Mpc]$')
    for irow in [irow_dens_fit]:
        axs[irow].set_yscale("log")
    
    for irow in range(nrows):
        for jrow in range(ncols):
            axs[irow, jrow].tick_params(axis='both', which='both', direction='in', top=True, right=True, length=6)

    plt.tight_layout()

    path = '/cosma7/data/dp004/dc-wang4/VVV_plots/WarmHalo/HaloDensity/'+\
    'HaloDensity_{}.png'.format(label)
    print('save image: ', path)
    fig.savefig(path, bbox_inches = 'tight', dpi=300)
    plt.close(fig)



    # r_ell_list = dist_bin
    # for i, selection in enumerate([semi_selection+(~semi_selection), semi_selection, ava_selection]):
    #     par_ava_elli_selection = dist_halo[selection] < np.max(dist_bin)*100

    #     par_x_ava = par_x_inhalo[par_ava_elli_selection]
    #     par_y_ava = par_y_inhalo[par_ava_elli_selection]
    #     par_z_ava = par_z_inhalo[par_ava_elli_selection]
    #     par_dist_ava = dist_halo[par_ava_elli_selection]
        
    #     for j in range(len(r_ell_list)):
    #         r_ell = r_ell_list[j]/halo_r
    #         comp, eigenvalue, eli_selection, x_eli, y_eli, z_eli, coor_transformed = \
    #         halo_ellipse(halo_pos, par_x_ava, par_y_ava, par_z_ava, halo_r, r_ell)

    #         label_r = label + 'r{}s{}'.format(j, i)
    #         Ell_Comp[label_r] = [eigenvalue, comp, eli_selection, coor_transformed]
            
    #         if j == len(r_ell_list)-1:
    #             Ell_LargestR[label] = np.max(par_dist_ava[eli_selection])








    


Halo_data = []

Halo_data.append(Chosen_Peak)
Halo_data.append(Halo_CausticsR)
Halo_data.append(LocalRho_Halo)
Halo_data.append(LocalRho_Semi)
Halo_data.append(RhoFitParam_Halo)
Halo_data.append(RhoFitParam_Semi)
Halo_data.append(RhoFitChi2_Halo)
Halo_data.append(RhoFitChi2_Semi)
Halo_data.append(Spear_Coef)
Halo_data.append(Redshift)
Halo_data.append(HaloMass_Data)
Halo_data.append(LocalRho_Various)
Halo_data.append(Peak_Density)
Halo_data.append(Peak_Dist_toCentre)
Halo_data.append(Selections)

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Final_Data/HaloDensity_{}'.format(label)
print('save data: ', path)
np.save(path, Halo_data)

print('Done all: '+str((time.clock()-t1)/60)+' min')
