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

k_max = ks_l[max_delta]

print('k_max: ', k_max)

Hconst = 0.677

SubHalo_Mass = {}
Halo_FirSub = {}
Halo_NSub = {}
Subhalo_FirstProg = {}
Halo_Index = {}

for T_min_temp in range(0, T_min+1):
    label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min_temp, T_min_temp)

    path = '/cosma7/data/dp004/dc-wang4/VVV_data/Ref_Data/AllSubRef_{}.npy'.format(label)
    try:
        AllHalo_data = np.load(path, allow_pickle = True)
        
        Halo_Index.update(AllHalo_data[11])
        print(path)
    except:
        pass


for T in range(0, T_max+1):
    
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

        # Par_Coor[label] = File[0]/Hconst
        # Par_Mass[label] = File[1]*10**10/Hconst
        # Par_ID[label] = File[2]
        # Par_Vel[label] = File[3]/Hconst
        # Halo_Mass_200[label] = File[4]*10**10/Hconst
        SubHalo_Mass[label] = File[5]*10**10/Hconst
        Halo_FirSub[label] = File[6]
        Halo_NSub[label] = File[7]
        # SubHalo_GrN[label] = File[8]
        # Halo_Offset[label] = File[9]
        # Halo_ParN[label] = File[10]
        # SubHalo_Offset[label] = File[11]
        # SubHalo_ParN[label] = File[12]
        # Halo_Mass[label] = File[13]*10**10/Hconst
        # Halo_Coor[label] = File[14]/Hconst
        # Subhalo_Coor[label] = File[15]/Hconst
        # Halo_R_200[label] = File[16]/Hconst
        # Redshift[label] = File[17]
        # Box_size[label] = File[18]/Hconst
        # Halo_Vel[label] = File[19]/Hconst
        # Subhalo_Spin[label] = File[20]
        if T > 0:
            Subhalo_FirstProg[label] = Read_merger_tree(L, T, 16, Path)
    except:
        T_max = T-1
        break
    

# print('finished reading', T_min, T_max)
if T_max < T_min:
    print('no file found')
    exit()

print('finished reading', T_min, T_max)

Subhalo_ProgMass = {}
for T in range(T_min, T_max+1):
    label = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T)

    try:
        halo_index = int(Halo_Index[label])
    except:
        print('No progenitor at ', label)
        continue

    halo_firstsub = Halo_FirSub[label][halo_index]
    halo_nsub = Halo_NSub[label][halo_index]

    subhalo_most_mass = np.zeros(halo_nsub)

    for isub in range(halo_firstsub, halo_firstsub + halo_nsub):
        sub_mass = SubHalo_Mass[label][isub]
        subhalo_most_mass[isub - halo_firstsub] = sub_mass

        if sub_mass < 0:
            continue

        subhalo_mass_prog = []
        sub_temp = int(isub)

        for T_prog in range(T, 0, -1):
            label_prog   = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_prog)
            label_prog_1 = 'L{}HaloID{}k{}Nh{}S{}'.format(L, HaloID, kfs, Nh, T_prog - 1)

            try:
                subhalo_firstprog = int(Subhalo_FirstProg[label_prog][sub_temp])

                if subhalo_firstprog < 0:
                    break

                subhalo_mass_prog.append(SubHalo_Mass[label_prog_1][subhalo_firstprog])
                sub_temp = subhalo_firstprog

            except:
                break

        if len(subhalo_mass_prog) > 0:
            subhalo_most_mass[isub - halo_firstsub] = np.max(np.append(sub_mass, subhalo_mass_prog))
        else:
            subhalo_most_mass[isub - halo_firstsub] = sub_mass

    Subhalo_ProgMass[label] = subhalo_most_mass


SubHalo_Data = []
SubHalo_Data.append(Subhalo_ProgMass)

label = 'L{}HaloID{}k{}Nh{}Tmin{}Tmax{}'.format(L, HaloID, kfs, Nh, T_min, T_max)
path = '/cosma7/data/dp004/dc-wang4/VVV_data/Final_Data/_{}'.format(label)
print('save data: ', path)
np.save(path, SubHalo_Data)

print('Done all: '+str((time.clock()-t1)/60)+' min')

                
            


