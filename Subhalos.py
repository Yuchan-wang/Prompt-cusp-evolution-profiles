import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import warnings; warnings.simplefilter('ignore')
import numpy.linalg as lina
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
import pandas as pd
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

Peaks_Properties = {}
# CoorAllHaloPar0 = {}
# LinRhoAllHaloPar = {}
Peak_IDs = {}
Semi_Selection = {}

print('finish initializing')

# for T_min_temp in range(0, T_max+1):
    # for T_max_temp in range(200, T_min_temp, -1):
label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_min)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Peak_Data/AllPeak_{}.npy'.format(label)
try:
    Peak_data = np.load(path, allow_pickle = True)
    
    Peaks_Properties.update(Peak_data[0])
    # CoorAllHaloPar0.update(Peak_data[1])
    # LinRhoAllHaloPar.update(Peak_data[2])
    Peak_IDs.update(Peak_data[3])
    Semi_Selection.update(Peak_data[4])
    print(path)
except:
    print('No peak data!')
    exit()

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
    print('No peak identifier data!')
    pass

# Indexer_Halo = {}
# Halo_Dist = {}
# Dist_Axis = {}
# Halo_Density = {}
# Par_SubIndex = {}
Subhalo_FirstProg = {}
# SubHalo_Mass = {}
# Halo_Mass = {}
Halo_Index = {}

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

nrows = 3
ncols = 4
NbinH = 70
Nrand = int(1e3)
# dens_factorH = 100
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
# par_separation = func_par_separation(Par_x0)
high_selection = Par_m0 == Par_m0[0]
Par_x0_m = Par_x0[high_selection]
Par_y0_m = Par_y0[high_selection]
Par_z0_m = Par_z0[high_selection]
range_x = np.max(Par_x0_m) - np.min(Par_x0_m)
range_y = np.max(Par_y0_m) - np.min(Par_y0_m)
range_z = np.max(Par_z0_m) - np.min(Par_z0_m)
par_separation = np.mean([range_x, range_y, range_z])/(Nh*20)


Halo_PeakCoor = {}
Halo_PeakID = {}
Halo_Peak_Rs = {}
Halo_Mr_Peak = {}
Halo_PeakLin_bin = {}
Halo_meanlinrho_Peak = {}
# Halo_EllParams = {}
# Halo_PeakEllR = {}
Halo_ElliFit = {}
Halo_PeakQua = {}
Halo_PeakElli = {}
Halo_PeakProl = {}
Halo_PeakPar = {}
Sub_Info = {}
Par_Peak = {}
# SubHalo_Prog_T = {}
# SubHalo_Prog_Index = {}
# SubHalo_MassiveProg = {}

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
        SubHalo_Mass[label] = File[5]*10**10/Hconst
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
    except:
        T_max = T-1
        break
    

# print('finished reading', T_min, T_max)
if T_max < T_min:
    print('no file found')
    exit()

print('finished reading', T_min, T_max)

mean_density = np.sum(Par_m0)/Box_size[label]**3
halo_mlim = SpuriousHalo_Mlim(mean_density, par_separation, k_max*Hconst)
print('{:e}'.format(halo_mlim))
label_L = 'L{}HaloID{}k{}Nh{}'.format(L, HaloID, kfs, Nh)
Halo_Mlim = {}
Halo_Mlim[label_L] = halo_mlim

halo_index = 0
r_high_lim = halor
N_nfw = 51

# SubPeak_Index = {}
SubPeak_Score = {}
SubHalo_ParCoor0 = {}
SubHalo_ParCoor = {}
SubHalo_LagCA = {}
SubHalo_MassList = {}
AllLinRhoHalo = {}
Paired_PeakSub = {}
Halo_Sub_List = {}

rfactor = 12/7
for T in range(T_min, T_max+1):    
    plt.clf()
    fig, axs = plt.subplots(figsize = (ncols*2+4, nrows*2+4),  #sharex = 'row', sharey = 'row',
                            nrows = nrows, ncols = ncols, 
                            dpi = 200)
    axs = axs.flat
    
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)

    if len(Halo_Mass_200[label]) == 0:
        print('No Halo at ', label)
        continue

    try:
        halo_index = int(Halo_Index[label])
    except:
        print('No progenitor at ', label)
        continue
    
    try:
        peak_ids = Peak_IDs[label]
        # npeak = Peaks_Properties[label]
        npeak = int(np.max(NPeak_Test[label][1]))
    except:
        print('Peak error!')
        continue
    if npeak == 0:
        print('No peak found!')
        continue
    index_peak = pd.Index(peak_ids)
    
    rfs = np.pi/kfs
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
        # halo_index = int(halo_prog_index[T])
        # halo_index = int(0)
        halo_x = Halo_Coor[label][halo_index, 0]
        halo_y = Halo_Coor[label][halo_index, 1]
        halo_z = Halo_Coor[label][halo_index, 2]
        halo_r = Halo_R_200[label][halo_index]
        halo_firstsub = Halo_FirSub[label][halo_index]
        halo_pos = [halo_x, halo_y, halo_z]
        halo_nsub = Halo_NSub[label][halo_index]


        if halo_nsub == 0:
            print ('no subhalo in:', label)
            continue
        print('{} subhalos in halo'.format(halo_nsub))

        sub_index_high = halo_firstsub + halo_nsub
        halo_sub_list = range(halo_firstsub, sub_index_high)
        Halo_Sub_List[label] = halo_sub_list

        halo_np = Halo_ParN[label][halo_index]
        halo_off = Halo_Offset[label][halo_index]
        
        # par_selection = func_SphSelect(Par_x, Par_y, Par_z, halo_r*halor,
        #                             halo_x, halo_y, halo_z)
    
        # par_x_inhalo = Par_x[par_selection]
        # par_y_inhalo = Par_y[par_selection]
        # par_z_inhalo = Par_z[par_selection]
        # par_m_inhalo = Par_m[par_selection]
        # par_id_inhalo = Par_id[par_selection]
        par_x_inhalo = Par_x[int(halo_off):int(halo_off+halo_np)]
        par_y_inhalo = Par_y[int(halo_off):int(halo_off+halo_np)]
        par_z_inhalo = Par_z[int(halo_off):int(halo_off+halo_np)]
        # par_m_inhalo = Par_m[int(halo_off):int(halo_off+halo_np)]
        par_id_inhalo = Par_id[int(halo_off):int(halo_off+halo_np)]
        
        par_selection = func_BoxSelect(Par_x, Par_y, Par_z,
                                    np.min(par_x_inhalo), np.max(par_x_inhalo), 
                                    np.min(par_y_inhalo), np.max(par_y_inhalo), 
                                    np.min(par_z_inhalo), np.max(par_z_inhalo))
        
        par_x_inbox = Par_x[par_selection]
        par_y_inbox = Par_y[par_selection]

        dens_map, xedges, yedges = np.histogram2d(par_x_inbox, par_y_inbox, bins=[NbinH, NbinH]) 
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        dens_factorH = np.mean(dens_map)
        axs[0].imshow(np.transpose(np.log10(dens_map+dens_factorH)), 
                                extent=extent, origin='lower',
                                cmap=cm.twilight_shifted, interpolation="spline16")
        
        par_inhalo_selection = index_parent.get_indexer(par_id_inhalo)
        par_inhalo_selection = par_inhalo_selection[par_inhalo_selection >= 0]

        par_x_inhalo0 = Par_x0[par_inhalo_selection]
        par_y_inhalo0 = Par_y0[par_inhalo_selection]
        par_z_inhalo0 = Par_z0[par_inhalo_selection]
        halo_mean_centre0 = [np.mean(par_x_inhalo0), 
                            np.mean(par_y_inhalo0), np.mean(par_z_inhalo0)]
        par_linrho_inhalo = Par_linrho[par_inhalo_selection]
        AllLinRhoHalo[label] = par_linrho_inhalo
        print('max linrho: ', np.max(AllLinRhoHalo[label]))
        print('linrho len: ', len(AllLinRhoHalo[label]))

        for iax in [1, 2, 3]:
            func_scatter(axs[iax], par_x_inhalo0, par_y_inhalo0, Nrand, 0.5, 0.5, c = 'grey')
            # func_scatter(axs[2], par_x_inhalo0, par_linrho_inhalo, Nrand, 0.5, 0.5, c = 'grey')
        
        # ntry = np.min([halo_nsub, ])
        ntry = halo_nsub
        normal = colors.Normalize(vmin = 0, vmax = ntry)
        
        
        normal_peak = colors.Normalize(vmin = 0, vmax = npeak)

        all_pair_scores = []
        for ipeak in range(0, npeak):
            label_peak = '{}P{}'.format(label, ipeak)
            peak_overlap = []
            halo_index = []
            
            try:
                # [centre_maxlinrho_temp, peak_maxlinrho_temp, 
                # x_axis, linrho_near_maxlinrho_all, 
                # quad_coef, log_coefs, ell_fit_coef, 
                # peak_r, ellipticity, prolateness, 
                # id_near_maxlinrho_temp,
                # x_near_maxlinrho_all, 
                # y_near_maxlinrho_all, 
                # z_near_maxlinrho_all] = Peaks_Properties[label_peak]
                [centre_maxlinrho_temp, peak_maxlinrho_temp, 
                x_axis, linrho_near_maxlinrho_all, 
                quad_coef, log_coefs, ell_fit_coef, 
                peak_r, ellipticity, prolateness, 
                eigenvals, eigenvecs, ell_fit_matrix,
                # ca,eigenvals
                id_near_maxlinrho_temp,
                x_near_maxlinrho_all, 
                y_near_maxlinrho_temp, 
                z_near_maxlinrho_temp] = Peaks_Properties[label_peak]
                semi_selection = Semi_Selection[label_peak]

                print('looking at peak: ', label_peak)
                print('Peak info:', peak_maxlinrho_temp, centre_maxlinrho_temp)

            except:
                print('peak unavailable', label_peak)
                continue
                
            x_near_maxlinrho_temp = par_x_inhalo0[semi_selection]
            y_near_maxlinrho_temp = par_y_inhalo0[semi_selection]
            z_near_maxlinrho_temp = par_z_inhalo0[semi_selection]
            linrho_near_maxlinrho_temp = par_linrho_inhalo[semi_selection]

            if len(x_near_maxlinrho_temp) > Nrand:
                index_array = np.array(range(len(x_near_maxlinrho_temp)))
                rand_index = np.random.choice(index_array, Nrand, replace=False)

                x_near_maxlinrho_temp_reduced = x_near_maxlinrho_temp[rand_index]
                y_near_maxlinrho_temp_reduced = y_near_maxlinrho_temp[rand_index]
                z_near_maxlinrho_temp_reduced = z_near_maxlinrho_temp[rand_index]
                linrho_near_maxlinrho_temp_reduced = linrho_near_maxlinrho_temp[rand_index]
            else:
                x_near_maxlinrho_temp_reduced = x_near_maxlinrho_temp
                y_near_maxlinrho_temp_reduced = y_near_maxlinrho_temp
                z_near_maxlinrho_temp_reduced = z_near_maxlinrho_temp
                # par_linrho_insub_temp_reduced = par_linrho_insub_temp
                linrho_near_maxlinrho_temp_reduced = linrho_near_maxlinrho_temp

            for isub in range(ntry):
                
                label_sub = 'L{}HaloID{}k{}Nh{}S{}SH{}'.format(L, HaloID, kfs, Nh, T, isub)
                sub_index = int(halo_firstsub + isub)
                print('looking at {}th subhalo'.format(sub_index))
                par_insub_off = SubHalo_Offset[label][int(sub_index)]
                par_insub_np = SubHalo_ParN[label][int(sub_index)]
                subhalo_mass = SubHalo_Mass[label][int(sub_index)]

                print('{} particles in subhalo'.format(par_insub_np))
                # m_sub = SubHalo_M

                
                
                par_id_insub = Par_id[int(par_insub_off):int(par_insub_off + par_insub_np)]
                par_x_insub_temp = Par_x[int(par_insub_off):int(par_insub_off + par_insub_np)]
                par_y_insub_temp = Par_y[int(par_insub_off):int(par_insub_off + par_insub_np)]
                par_z_insub_temp = Par_z[int(par_insub_off):int(par_insub_off + par_insub_np)]

                par_off_low = int(par_insub_off - halo_off)
                par_off_high = int(par_insub_off - halo_off + par_insub_np)
                # par_insub_selection = index_parent.get_indexer(par_id_insub)
                # par_insub_selection = par_insub_selection[par_insub_selection >= 0]
                par_x_insub0_temp = par_x_inhalo0[par_off_low:par_off_high]
                par_y_insub0_temp = par_y_inhalo0[par_off_low:par_off_high]
                par_z_insub0_temp = par_z_inhalo0[par_off_low:par_off_high]
                par_linrho_insub_temp = par_linrho_inhalo[par_off_low:par_off_high]

                if ipeak == 0:
                    func_scatter(axs[1], par_x_insub0_temp, par_y_insub0_temp, 
                            Nrand, 0.5, 0.5, c = cm.rainbow(normal(isub)))
                    
                    x_c = par_x_insub0_temp - np.mean(par_x_insub0_temp)
                    y_c = par_y_insub0_temp - np.mean(par_y_insub0_temp)
                    z_c = par_z_insub0_temp - np.mean(par_z_insub0_temp)
                    a_axis = b_axis = c_axis = 1.0
                    eigenvec_par = np.identity(3)

                    for _ in range(10):
                        inertia = reduced_inertia(x_c, y_c, z_c, a_axis, b_axis, c_axis, eigenvec_par)
                        eigenval_par, eigenvec_par = lina.eigh(inertia)
                        axes = np.sqrt(np.abs(eigenval_par))
                        c_axis = axes[0]
                        b_axis = axes[1]
                        a_axis = axes[2]


                    # inertia = reduced_inertia(par_x_insub0_temp, par_y_insub0_temp, par_z_insub0_temp)
                    eigenval_par, eigenvec_par =  lina.eigh(inertia)
                    ca = np.min(eigenval_par)/np.max(eigenval_par)

                    # func_scatter(axs[2], par_x_insub0_temp, par_y_insub0_temp, 
                    #             Nrand, 0.5, 0.5, c = cm.rainbow((normal(isub))))  
                    

                    SubHalo_LagCA[label_sub] = ca
                    SubHalo_MassList[label_sub] = subhalo_mass
                    SubHalo_ParCoor0[label_sub] = [par_x_insub0_temp, par_y_insub0_temp, par_z_insub0_temp]
                    SubHalo_ParCoor[label_sub] = [par_x_insub_temp, par_y_insub_temp, par_z_insub_temp]
                
                # dist_sub = distance(par_x_insub0_temp, par_y_insub0_temp, par_z_insub0_temp, 
                #                     halo_mean_centre0)
                # dist_sub_mean = np.mean(dist_sub)
                # dist_sub_std = np.std(dist_sub)

                # halo_box0_min = [np.min(par_x_insub0_temp), 
                #                 np.min(par_y_insub0_temp), 
                #                 np.min(par_z_insub0_temp)]
                # halo_box0_max = [np.max(par_x_insub0_temp), 
                #                 np.max(par_y_insub0_temp), 
                #                 np.max(par_z_insub0_temp)]
                # halo_box_min = [np.min(par_x_insub_temp), 
                #                 np.min(par_y_insub_temp), 
                #                 np.min(par_z_insub_temp)]
                # halo_box_max = [np.max(par_x_insub_temp), 
                #                 np.max(par_y_insub_temp), 
                #                 np.max(par_z_insub_temp)]

                
            
                # isinsub0 = IsInBox(centre_maxlinrho_temp, halo_box0_min, halo_box0_max)
                # isinsub = IsInBox(peak_maxlinrho_temp, halo_box_min, halo_box_max)
                # if (isinsub0 == False) or (isinsub == False):
                #     # continue
                #     pass
                

                # print('before, linrho: ', np.min(par_linrho_insub_temp_reduced), 
                #     np.max(par_linrho_insub_temp_reduced))
                
                # print('after, linrho: ', np.min(par_linrho_insub_temp_reduced), 
                #     np.max(par_linrho_insub_temp_reduced))
                # else:
                if len(par_x_insub0_temp) > Nrand:
                    index_array = np.array(range(len(par_x_insub0_temp)))
                    rand_index = np.random.choice(index_array, Nrand, replace=False)
                    par_x_insub0_temp_redeuced = par_x_insub0_temp[rand_index]
                    par_y_insub0_temp_redeuced = par_y_insub0_temp[rand_index]
                    par_z_insub0_temp_redeuced = par_z_insub0_temp[rand_index]
                    par_linrho_insub_temp_reduced = par_linrho_insub_temp[rand_index]
                else:
                    par_x_insub0_temp_redeuced = par_x_insub0_temp
                    par_y_insub0_temp_redeuced = par_y_insub0_temp
                    par_z_insub0_temp_redeuced = par_z_insub0_temp
                    par_linrho_insub_temp_reduced = par_linrho_insub_temp

                par_linrho_insub_temp_reduced += 1
                linrho_near_maxlinrho_temp_reduced += 1
                print('length:', len(x_near_maxlinrho_temp_reduced), 
                len(par_x_insub0_temp_redeuced), len(par_linrho_insub_temp_reduced))
                score = Score_Overlap(x_near_maxlinrho_temp_reduced,
                                    y_near_maxlinrho_temp_reduced, 
                                    z_near_maxlinrho_temp_reduced, 
                                    par_x_insub0_temp_redeuced, 
                                    par_y_insub0_temp_redeuced, 
                                    par_z_insub0_temp_redeuced, 
                                    linrho_near_maxlinrho_temp_reduced,
                                    par_linrho_insub_temp_reduced,
                                    par_separation)
                label_peaksub = '{}P{}SH{}'.format(label, ipeak, isub)
                
                print(label_peaksub, score)
                if score == 0:
                    continue
                score = np.log10(score)
                    
                peak_overlap.append(score)
                halo_index.append(isub)
                
                SubPeak_Score[label_peaksub] = score
                all_pair_scores.append((score, ipeak, isub))
                # print('label_peaksub', score)

            # SubPeak_Index[label_peak] = halo_index
            # SubPeak_Score[label_peak] = peak_overlap

            peak_overlap = np.array(peak_overlap)
            halo_index = np.array(halo_index)
            print('Peak Score', len(halo_index), len(peak_overlap), label_sub)
            
            if len(halo_index) == 0:
                npeak_insub = 0
                color_peak = 'green'


            else:

                all_pair_scores.sort(key=lambda x: x[0], reverse=True)

                # axs[7].scatter(peak_index, peak_overlap, color = cm.Greys(normal(isub)), 
                #                , marker = 'x', s = 100)
                # axs[2].scatter(peak_index, peak_overlap,
                #             color = cm.Greys(normal(isub)), marker = '.', s = 20)
                peak_available_selection = peak_overlap > -1000
                # npeak_insub = np.sum(peak_available_selection)
                nsub_inpeak = np.sum(peak_available_selection)
                peak_overlap_available = peak_overlap[peak_available_selection]
                halo_index_available = halo_index[peak_available_selection]


                matched_peaks = set()
                matched_subs = set()
                Paired_PeakSub = {}


                if nsub_inpeak == 0:
                    color_peak = 'green'
                    
                else:
                    for score, p_idx, s_idx in all_pair_scores:
                        if p_idx not in matched_peaks and s_idx not in matched_subs:
                            matched_peaks.add(p_idx)
                            matched_subs.add(s_idx)
                            label_sub_temp = '{}SH{}'.format(label, s_idx)
                            Paired_PeakSub[label_sub_temp] = p_idx
                            label_peak_temp = '{}P{}'.format(label, p_idx)
                            Paired_PeakSub[label_peak_temp] = s_idx
                    
                    # peak_insub_index = np.argmax(peak_overlap_available)
                    # score_insub = peak_overlap_available[peak_insub_index]
                    # isub_inpeak = halo_index_available[peak_insub_index]

                    # Paired_PeakSub[label_peak] = isub_inpeak
                    # label_sub_temp = '{}SH{}'.format(label, isub_inpeak)
                    # try:
                    #     existing_pair = Paired_PeakSub[label_sub_temp]
                    #     existing_label = '{}P{}SH{}'.format(label, existing_pair, isub_inpeak)
                    #     new_label = '{}P{}SH{}'.format(label, ipeak, isub_inpeak)
                    #     if SubPeak_Score[existing_label] > SubPeak_Score[new_label]:
                    #         pass
                    #     else:
                    #         Paired_PeakSub[label_sub_temp] = ipeak
                    # except:
                    #     Paired_PeakSub[label_sub_temp] = ipeak
                    # print('Found {} and {} as a pair'.format(label_peak, isub_inpeak))
                    color_peak = 'red'
                    # 
                    # color_sub = 'red'
                    # axs[8].scatter(ipeak, score_insub,
                    #         color = cm.rainbow(normal_peak(isub)), marker = 'x', s = 100)

            # print('{} particles as peak'.format(npeak_insub))

            axs[3].scatter(centre_maxlinrho_temp[0], centre_maxlinrho_temp[1], 
                               color = color_peak,
                               marker = 'x', s = 200)
            
            for isub in range(0, ntry):
                label_sub = 'L{}HaloID{}k{}Nh{}S{}SH{}'.format(L, HaloID, kfs, Nh, T, isub)
                score_list = []
                score_max = 0
                try:
                    subhalo_mass = SubHalo_MassList[label_sub]
                    [par_x_insub0_temp, par_y_insub0_temp, par_z_insub0_temp] = SubHalo_ParCoor0[label_sub]
                    ca = SubHalo_LagCA[label_sub]
                except:

                    continue

                try:
                    ipeak_paired = Paired_PeakSub[label_sub]
                    color_sub = 'red'
                    
                except:
                    color_sub = 'green'

                for ipeak in range(npeak):
                    label_peaksub = '{}P{}SH{}'.format(label, ipeak, isub)
                    try:
                        score_list.append(SubPeak_Score[label_peaksub])
                    except:
                        continue
                    
                    score_list = np.array(score_list)
                    score_max = np.max(score_list)
                
                func_scatter(axs[2], par_x_insub0_temp, par_y_insub0_temp, 
                                Nrand, 0.5, 0.5, c = color_sub)
                axs[4].scatter(subhalo_mass, ca, 
                            color = color_sub, s = 100, marker = 'x', alpha = 0.5)
                axs[5].scatter(subhalo_mass, score_max, 
                            color = color_sub, s = 100, marker = 'x', alpha = 0.5)
                axs[6].scatter(score_max, ca, 
                            color = color_sub, s = 100, marker = 'x', alpha = 0.5)

            # prog_T_list = []
            # prog_index_list = []
            # prog_mass_list = []
            # prog_index = sub_index
            # prog_max_mass = subhalo_mass
            # for T_prog in range(T, 0, -1):
            #     prog_index_list.append(prog_index)
            #     prog_T_list.append(T_prog)

                
            #     label_prog = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_prog)
            #     try:
            #         prog_mass = SubHalo_Mass[label_prog][prog_index]
                    
            #         # if prog_mass > 3*prog_max_mass:
            #         #     break
            #         prog_mass_list.append(prog_mass)

            #         if prog_mass > prog_max_mass:
            #             prog_max_mass = prog_mass

            #         prog_index = Subhalo_FirstProg[label_prog][prog_index]
            #         print(label_sub, label_prog, prog_index)
            #     except:
            #         break
                
            #     if prog_index == -1:
            #         break
            # SubHalo_Prog_T[label_sub] = prog_T_list
            # SubHalo_Prog_Index[label_sub] = prog_index_list

            # prog_mass_list = np.array(prog_mass_list)
            # SubHalo_MassiveProg[label_sub] = [prog_T_list, prog_index_list, prog_mass_list, np.max(prog_mass_list)]

            
            # fractions = []
            # z_list = []
            # normal_T = colors.Normalize(vmin = T, vmax = T_max)
            # for T_ref in range(T, T_max+1):
            #     label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
            #     index_halo = Indexer_Halo[label]
            #     issubhaloin = index_halo.get_indexer(par_id_insub)

            #     npar_in = np.sum(issubhaloin >= 0)
            #     issubhaloin = issubhaloin[issubhaloin >= 0]

            #     fraction_in = npar_in/len(par_id_insub)
            #     z_list.append(Redshift[label])
            #     fractions.append(fraction_in)
            #     print(label, '{} found for {}th subhalo'.format(fraction_in, isub))

            #     if npar_in > 0:
            #         dist_halo = Halo_Dist[label]
            #         dist_halo_fromsub = dist_halo[issubhaloin]
            #         mean_dist_fromsub = np.mean(dist_halo_fromsub)
            #         std_dist_fromsub = np.std(dist_halo_fromsub)

            #         [dist_bin, d_dist, dist_bin_mid] = Dist_Axis[label]

            #         subdist_bins, subbin_edges = np.histogram(dist_halo_fromsub, 
            #                                             bins = dist_bin)
            #         density_sub = subdist_bins/d_dist/dist_bin_mid**2/4/np.pi*par_m_inhalo[0]
            #         if color_sub == 'red':
            #             Density_AllSub_Real[label] += density_sub
            #         else:
            #             Density_AllSub_Spurious[label] += density_sub


            #         dist_bin_mid = dist_bin_mid[density_sub > 0]
            #         density_sub = density_sub[density_sub > 0]
            #         axs[8].plot(dist_bin_mid, 
            #                     density_sub*dist_bin_mid**rfactor, 
            #                     color = cm.rainbow(normal_T(T_ref)), 
            #                     alpha = 0.6, lw = 1)
            #         axs[9].plot(dist_bin_mid, 
            #                     density_sub*dist_bin_mid**rfactor, 
            #                     color = color_sub, 
            #                     alpha = 0.6, lw = 1)
                    
                        # axs[10].plot(dist_bin_mid, 
                        #         density_sub*dist_bin_mid**rfactor, 
                        #         color = cm.rainbow(score_insub), 
                        #         alpha = 0.6, lw = 1)            
                    # bin_mid = bin_edges[1:] + bin_edges[:-1]
                    # axs[7].errorbar(dist_sub_mean, mean_dist_fromsub, 
                    #                 xerr = dist_sub_std, yerr = std_dist_fromsub, 
                    #                 color = cm.rainbow(normal_T(T_ref)), 
                    #                 alpha = 0.5)
                    # axs[8].errorbar(dist_sub_mean, mean_dist_fromsub, 
                    #                 xerr = dist_sub_std, yerr = std_dist_fromsub, 
                    #                 color = color_sub, 
                    #                 alpha = 0.5)
                    # axs[1, isub+1].plot(bin_mid, dist_bins/bin_mid**2, color = cm.rainbow(normal_T(T_ref)))
            
            
            # axs[5].plot(z_list, fractions, color = cm.rainbow(normal(isub)), alpha = 0.5)
            # axs[6].plot(z_list, fractions, color = color_sub, alpha = 0.5)
            
            # Sub_Info[label] = []

        


        # for T_ref in range(T+1, T_max+1):
        #     label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_ref)
        #     [dist_bin_mid, density_halo] = Halo_Density[label]
        #     for ipanel in [8, 9, 10]:
        #         axs[ipanel].plot(dist_bin_mid, 
        #                     density_halo*dist_bin_mid**rfactor, 
        #                     color = cm.rainbow(normal_T(T_ref)), linestyle = 'dotted',
        #                     alpha = 0.6, lw = 1)

        #     density_real = Density_AllSub_Real[label]
        #     density_spurious = Density_AllSub_Spurious[label]
        #     axs[10].plot(dist_bin_mid[density_real>0], 
        #                 density_real[density_real>0]*dist_bin_mid[density_real>0]**rfactor, 
        #                 color = cm.hot(normal_T(T_ref)), linestyle = 'solid',
        #                 alpha = 0.6, lw = 1)
        #     axs[10].plot(dist_bin_mid[density_spurious>0], 
        #                 density_spurious[density_spurious>0]*dist_bin_mid[density_spurious>0]**rfactor, 
        #                 color = cm.cool(normal_T(T_ref)), linestyle = 'solid',
        #                 alpha = 0.6, lw = 1)
        
        # for irow in range(nrows):
        #     for icol in range(ncols):
        for ipanel in range(nrows*ncols):
                axs[ipanel].tick_params(axis='both', which='both', direction='in', top=True, right=True, length=6)
        
        axs[4].axvline(x = halo_mlim, linestyle = 'dashed', color = 'Grey')
        axs[5].axvline(x = halo_mlim, linestyle = 'dashed', color = 'Grey')

        axs[4].axhline(y = 0.2, linestyle = 'dashed', color = 'Grey')
        axs[6].axhline(y = 0.2, linestyle = 'dashed', color = 'Grey')

        axs[5].axhline(y = 0.2, linestyle = 'dashed', color = 'Grey')
        axs[6].axvline(x = 0.2, linestyle = 'dashed', color = 'Grey')

        axs[4].set_xscale("log")
        axs[5].set_xscale("log")
        # axs[6].set_xscale("log")
        for ipanel in [8, 9, 10]:
            axs[ipanel].set_xscale("log")
            axs[ipanel].set_yscale("log")

        label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)
        plt.tight_layout()
        path = '/cosma7/data/dp004/dc-wang4/VVV_plots/WarmHalo/SubPeak/'+\
        'AllSubHalo_{}.png'.format(label)
        print('save image: ', path)
        fig.savefig(path, bbox_inches = 'tight', dpi=300)
        plt.close(fig)    

AllSub_data = []
# AllSub_data.append(SubPeak_Index)
AllSub_data.append(SubPeak_Score)
AllSub_data.append(SubHalo_ParCoor0)
AllSub_data.append(SubHalo_ParCoor)
AllSub_data.append(SubHalo_LagCA)
AllSub_data.append(SubHalo_MassList)
AllSub_data.append(Halo_Mlim)
AllSub_data.append(Paired_PeakSub)
AllSub_data.append(Halo_Sub_List)

# AllSub_data.append(SubHalo_Prog_T)
# AllSub_data.append(SubHalo_Prog_Index)
# AllSub_data.append(SubHalo_MassiveProg)


# AllSub_data.append(AllLinRhoHalo)

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/SubPeak_Data/AllSubPeak_{}'.format(label)
print('save data: ', path)
np.save(path, AllSub_data)

print('Done all: '+str((time.clock()-t1)/60)+' min')