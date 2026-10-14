from mpi4py import MPI
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()


import sys

#load modules
print('loading modules...')
import numpy as np
import basic
import pickle
import curvedsky
#import cmb
import healpy as hp
import pymaster as nmt
import matplotlib.pyplot as plt

##########################################################################

import os,sys
import pylab as pl
import numpy as np
import lenspyx
import healpy as hp

sys.path.append('/sharefs/alicpt/users/chenwz/download/cmblensplus2/utils')
sys.path.append('/sharefs/alicpt/users/chenwz/Testarea/FORECAST/Foreground/run_with_foreground')

from lensingb_mine import *
from noise_making import *
from utils_mine import *
from phi_noise_making import *
from generate_map import *
from lenspyx.utils import camb_clfile
from lenspyx.utils_hp import synalm, almxfl, alm2cl
from generate_gaussian_fg import *

camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))    
sys.path.insert(0,camb_path)

import camb
from camb import model, initialpower ,correlations

from scipy.special import factorial

import pymaster as nm
import random




# 获取传递给脚本的参数
args = sys.argv[1:]

# 打印当前作业的编号
# 计算当前作业的文件编号列表
file_indices = comm_rank + int(args[0]) #+ 200

# 打印当前作业的文件编号列表
print(f"Job {comm_rank} will process files: {file_indices}")


# map输出路径 
#QUP_filename = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/QUP/QUP_{file_indices}.fits"

#TQU_N_filename = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK/sim/noise/LAT/280/QUP_N_{file_indices}.npy"




#FG_filename_27 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/27/IQU_FG_{file_indices}.fits"
#FG_filename_39 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/39/IQU_FG_{file_indices}.fits"
#FG_filename_93 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/93/IQU_FG_{file_indices}.fits"
#FG_filename_145 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/145/IQU_FG_{file_indices}.fits"
#FG_filename_225 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/225/IQU_FG_{file_indices}.fits"
#FG_filename_280 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/280/IQU_FG_{file_indices}.fits"



#print(f"Job {comm_rank} will write map to {TQU_N_filename2 }:")






#phi noise
#nlev_phi = [5.8,6.3]
#theta_ac_phi = [2.2,1.4]    #uK-arcmin
#phi_noise = gaussian_phi_n0(lmin_cut=2,lmax_cut=3*nside-1,rlmin_cut=[200],rlmax_cut=[4800],nside=nside,sigmas= nlev_phi ,thetas=theta_ac_phi,flag=5,seed = np.random.randint(99999999, size=1))


#TQU noise
# Ali-SAT  : 27GHz, 39GHz, 93GHz, 145GHz, 225GHz, 280GHz

nside = 2048
lmax = 3*nside-1 
lmax_len = lmax # desired lmax of the lensed field.
dlmax = lmax    # lmax of the unlensed fields is lmax + dlmax.  (some buffer is required for accurate lensing at lmax)
epsilon = 1e-6  # target accuracy of the output maps (execution time has a fairly weak dependence on this)
lmax_unl, mmax_unl = lmax_len + dlmax, lmax_len + dlmax
Tcmb = 2.7255e6
npix = hp.nside2npix(nside)

freqs = [27, 40, 90, 150, 225, 280]
nlev = [18.27, 15.51, 3.60, 4.57, 11.46, 22.49]
theta_ac = [97.34, 65.70, 16.20, 9.70, 10.60, 8.70]

N_red_T_s = [100, 39, 230, 1500, 17000, 31000]
N_red_T = [ 0 * value / (8760 * 3600) for value in N_red_T_s]   # ALI SAT DOESN'T USE T FIELD
l_knee_T = [1000, 1000, 1000, 1000, 1000, 1000] 
alpha_knee_T = [-3.5, -3.5, -3.5, -3.5, -3.5, -3.5, ]

N_red_P = nlev
l_knee_P = [30.00, 30.00, 50.00, 50.00, 70.00, 100.00] 
alpha_knee_P = [-2.40, -2.40, -2.50, -3.00, -3.00, -3.00]

ac2rad = np.pi/10800.
theta_rad = np.array(theta_ac)*ac2rad

mask_ali = hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/mask/ALI_LAT_add_SO_LAT_cross_PLANCK_{nside}_bi.fits')

for i in range(len(freqs)):
    print(f"Making noise map for {freqs[i]}GHz...")
    noise = obs_noise(nlev[i], nside,  seed=random.randint(0, 99999999), N_red_T=N_red_T[i], l_knee_T=l_knee_T[i], alpha_knee_T=alpha_knee_T[i], n_red_P=N_red_P[i], l_knee_P=l_knee_P[i], alpha_knee_P=alpha_knee_P[i], if_red=1)
    
    #noise = homo_noise(nlev[i], nside, seed=random.randint(0, 99999999), fwhm_f=None)

    #noise *= mask_ali
    
    TQU_N_filename2 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/Ali/SAT/{freqs[i]}/TQU_N_{file_indices}.npy"

    TQU_N_filename2_root = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/Ali/SAT/{freqs[i]}"

    if not os.path.exists(TQU_N_filename2_root):
        os.makedirs(TQU_N_filename2_root)


    noise_compressed = map_compression(nside,map=[noise[0],noise[1],noise[2]],dat=None,savepath=TQU_N_filename2,flag='f2d')
#noise *= hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK/sim/mask/AliCPT_LAT_cross_PLANCK_{nside}_bi.fits')
#hp.write_map(f"{QU_N_filename}",[noise[1],noise[2]],)


#TQU noise
# Ali-LAT  : 27GHz, 39GHz, 93GHz, 145GHz, 225GHz, 280GHz

nside = 2048
lmax = 3*nside-1 
lmax_len = lmax # desired lmax of the lensed field.
dlmax = lmax    # lmax of the unlensed fields is lmax + dlmax.  (some buffer is required for accurate lensing at lmax)
epsilon = 1e-6  # target accuracy of the output maps (execution time has a fairly weak dependence on this)
lmax_unl, mmax_unl = lmax_len + dlmax, lmax_len + dlmax
Tcmb = 2.7255e6
npix = hp.nside2npix(nside)

freqs = [27, 40, 90, 150, 225, 280]
nlev = [52, 27, 5.8, 6.3, 15, 37]
theta_ac = [7.4, 5.1, 2.2, 1.4, 1.0, 0.9]

N_red_T_s = [100, 39, 230, 1500, 17000, 31000]
N_red_T = [value / (8760 * 3600) for value in N_red_T_s]
l_knee_T = [1000, 1000, 1000, 1000, 1000, 1000] 
alpha_knee_T = [-3.5, -3.5, -3.5, -3.5, -3.5, -3.5, ]

N_red_P = nlev
l_knee_P = [700, 700, 700, 700, 700, 700, ] 
alpha_knee_P = [-1.4, -1.4, -1.4, -1.4, -1.4, -1.4, ]

ac2rad = np.pi/10800.
theta_rad = np.array(theta_ac)*ac2rad

mask_ali = hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/mask/ALI_LAT_add_SO_LAT_cross_PLANCK_{nside}_bi.fits')

for i in range(len(freqs)):
    print(f"Making noise map for {freqs[i]}GHz...")
    noise = obs_noise(nlev[i], nside,  seed=random.randint(0, 99999999), N_red_T=N_red_T[i], l_knee_T=l_knee_T[i], alpha_knee_T=alpha_knee_T[i], n_red_P=N_red_P[i], l_knee_P=l_knee_P[i], alpha_knee_P=alpha_knee_P[i], if_red=1)
    
    #noise = homo_noise(nlev[i], nside, seed=random.randint(0, 99999999), fwhm_f=None)

    #noise *= mask_ali
    
    TQU_N_filename2 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/Ali/LAT/{freqs[i]}/TQU_N_{file_indices}.npy"

    TQU_N_filename2_root = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/Ali/LAT/{freqs[i]}"

    if not os.path.exists(TQU_N_filename2_root):
        os.makedirs(TQU_N_filename2_root)


    noise_compressed = map_compression(nside,map=[noise[0],noise[1],noise[2]],dat=None,savepath=TQU_N_filename2,flag='f2d')
#noise *= hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK/sim/mask/AliCPT_LAT_cross_PLANCK_{nside}_bi.fits')
#hp.write_map(f"{QU_N_filename}",[noise[1],noise[2]],)





"""
#TQU noise
# Ali-LAT  : 27GHz, 39GHz, 93GHz, 145GHz, 225GHz, 280GHz

nside = 2048
lmax = 3*nside-1 
lmax_len = lmax # desired lmax of the lensed field.
dlmax = lmax    # lmax of the unlensed fields is lmax + dlmax.  (some buffer is required for accurate lensing at lmax)
epsilon = 1e-6  # target accuracy of the output maps (execution time has a fairly weak dependence on this)
lmax_unl, mmax_unl = lmax_len + dlmax, lmax_len + dlmax
Tcmb = 2.7255e6
npix = hp.nside2npix(nside)

freqs = [27, 39, 93, 145, 225, 280]
nlev = [52, 27, 5.8, 6.3, 15, 37]
theta_ac = [7.4, 5.1, 2.2, 1.4, 1.0, 0.9]

N_red_T_s = [100, 39, 230, 1500, 17000, 31000]
N_red_T = [value / (8760 * 3600) for value in N_red_T_s]
l_knee_T = [1000, 1000, 1000, 1000, 1000, 1000] 
alpha_knee_T = [-3.5, -3.5, -3.5, -3.5, -3.5, -3.5, ]

N_red_P = nlev
l_knee_P = [700, 700, 700, 700, 700, 700, ] 
alpha_knee_P = [-1.4, -1.4, -1.4, -1.4, -1.4, -1.4, ]

ac2rad = np.pi/10800.
theta_rad = np.array(theta_ac)*ac2rad

mask_ali = hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/mask/ALI_LAT_add_SO_LAT_cross_PLANCK_{nside}_bi.fits')

for i in range(len(freqs)):
    print(f"Making noise map for {freqs[i]}GHz...")
    noise = obs_noise(nlev[i], nside,  seed=random.randint(0, 99999999), N_red_T=N_red_T[i], l_knee_T=l_knee_T[i], alpha_knee_T=alpha_knee_T[i], n_red_P=N_red_P[i], l_knee_P=l_knee_P[i], alpha_knee_P=alpha_knee_P[i], if_red=1)
    
    #noise = homo_noise(nlev[i], nside, seed=random.randint(0, 99999999), fwhm_f=None)

    #noise *= mask_ali
    
    TQU_N_filename2 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/Ali/{freqs[i]}/TQU_N_{file_indices}.npy"

    TQU_N_filename2_root = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/Ali/{freqs[i]}"

    if not os.path.exists(TQU_N_filename2_root):
        os.makedirs(TQU_N_filename2_root)


    noise_compressed = map_compression(nside,map=[noise[0],noise[1],noise[2]],dat=None,savepath=TQU_N_filename2,flag='f2d')
#noise *= hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK/sim/mask/AliCPT_LAT_cross_PLANCK_{nside}_bi.fits')
#hp.write_map(f"{QU_N_filename}",[noise[1],noise[2]],)
"""

"""

# SO-LAT  : 27GHz, 39GHz, 93GHz, 145GHz, 225GHz, 280GHz

nside = 2048
lmax = 3*nside-1 
lmax_len = lmax # desired lmax of the lensed field.
dlmax = lmax    # lmax of the unlensed fields is lmax + dlmax.  (some buffer is required for accurate lensing at lmax)
epsilon = 1e-6  # target accuracy of the output maps (execution time has a fairly weak dependence on this)
lmax_unl, mmax_unl = lmax_len + dlmax, lmax_len + dlmax
Tcmb = 2.7255e6
npix = hp.nside2npix(nside)

freqs = [27, 39, 93, 145, 225, 280]
nlev = [52, 27, 5.8, 6.3, 15, 37]
theta_ac = [7.4, 5.1, 2.2, 1.4, 1.0, 0.9]

N_red_T_s = [100, 39, 230, 1500, 17000, 31000]
N_red_T = [value / (8760 * 3600) for value in N_red_T_s]
l_knee_T = [1000, 1000, 1000, 1000, 1000, 1000] 
alpha_knee_T = [-3.5, -3.5, -3.5, -3.5, -3.5, -3.5, ]

N_red_P = nlev
l_knee_P = [700, 700, 700, 700, 700, 700, ] 
alpha_knee_P = [-1.4, -1.4, -1.4, -1.4, -1.4, -1.4, ]

ac2rad = np.pi/10800.
theta_rad = np.array(theta_ac)*ac2rad

mask_so = hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/mask/SO_LAT_cross_PLANCK_{nside}_bi.fits')

for i in range(len(freqs)):
    print(f"Making noise map for {freqs[i]}GHz...")
    noise = obs_noise(nlev[i], nside,  seed=random.randint(0, 99999999), N_red_T=N_red_T[i], l_knee_T=l_knee_T[i], alpha_knee_T=alpha_knee_T[i], n_red_P=N_red_P[i], l_knee_P=l_knee_P[i], alpha_knee_P=alpha_knee_P[i], if_red=1)
    
    #noise = homo_noise(nlev[i], nside, seed=random.randint(0, 99999999), fwhm_f=None)
    #noise *= mask_so

    TQU_N_filename2 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/SO/{freqs[i]}/TQU_N_{file_indices}.npy"

    TQU_N_filename2_root = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise_red/SO/{freqs[i]}"
    if not os.path.exists(TQU_N_filename2_root):
        os.makedirs(TQU_N_filename2_root)


    noise_compressed = map_compression(nside,map=[noise[0],noise[1],noise[2]],dat=None,savepath=TQU_N_filename2,flag='f2d')









# LiteBIRD : 40GHz, 50GHz, 60GHz, 68GHz, 78GHz, 89GHz, 100GHz, 119GHz, 140GHz, 166GHz, 195GHz, 235GHz, 280GHz, 337GHz, 402GHz

nside = 128
lmax = 3*nside-1 
lmax_len = lmax # desired lmax of the lensed field.
dlmax = lmax    # lmax of the unlensed fields is lmax + dlmax.  (some buffer is required for accurate lensing at lmax)
epsilon = 1e-6  # target accuracy of the output maps (execution time has a fairly weak dependence on this)
lmax_unl, mmax_unl = lmax_len + dlmax, lmax_len + dlmax
Tcmb = 2.7255e6
npix = hp.nside2npix(nside)

freqs = [40, 50, 60, 68, 78, 89, 100, 119, 140, 166, 195, 235, 280, 337, 402]
theta_ac = [70.5, 58.5, 51.1, 41.6, 36.9, 33.0, 30.2, 26.3, 23.7, 28.9, 28.0, 24.7, 22.5, 20.9, 17.9,]
nlev = [26.5, 23.7, 15.1, 14.1, 11.0, 8.7, 7.3, 5.4, 5.1, 3.9, 5.0, 7.6, 9.8, 15.5, 33.6, ]

ac2rad = np.pi/10800.
theta_rad = np.array(theta_ac)*ac2rad

for i in range(len(freqs)):
    print(f"Making noise map for {freqs[i]}GHz...")
    #noise = obs_noise(nlev[i], nside,  seed=random.randint(0, 99999999), N_red_T=None, l_knee_T=None, alpha_knee_T=None, n_red_P=None, l_knee_P=None, alpha_knee_P=None, if_red=False)

    noise = homo_noise(nlev[i], nside, seed=random.randint(0, 99999999), fwhm_f=None)
    
    TQU_N_filename2 = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise/LiteBIRD/{freqs[i]}/TQU_N_{file_indices}.npy"

    TQU_N_filename2_root = f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_MAIN/sim/noise/LiteBIRD/{freqs[i]}"
    if not os.path.exists(TQU_N_filename2_root):
        os.makedirs(TQU_N_filename2_root)


    noise_compressed = map_compression(nside,map=[noise[0],noise[1],noise[2]],dat=None,savepath=TQU_N_filename2,flag='f2d')
#noise *= hp.read_map(f'/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK/sim/mask/AliCPT_LAT_cross_PLANCK_{nside}_bi.fits')
#hp.write_map(f"{QU_N_filename}",[noise[1],noise[2]],)


"""




















"""
#FG
nus= {27: 27*1e9, 39 : 39*1e9, 93 : 93*1e9, 145 : 145*1e9, 225 : 225*1e9, 280 : 280*1e9,}
fg_map,_ = gaussian_fg(2048, nus, file_indices + 9999, same=False, method=1)
hp.write_map(f"{FG_filename_27}",fg_map[27],)
hp.write_map(f"{FG_filename_39}",fg_map[39],)
hp.write_map(f"{FG_filename_93}",fg_map[93],)
hp.write_map(f"{FG_filename_145}",fg_map[145],)
hp.write_map(f"{FG_filename_225}",fg_map[225],)
hp.write_map(f"{FG_filename_280}",fg_map[280],)
"""

"""
#CMB 
Tlen_full,Qlen_full,Ulen_full,Pmap_full = making_maps_new(nside, nrms_f=None,fwhm_f=None,phi_map=None,pixwin=False)
hp.write_map(f"{QUP_filename}",[Qlen_full,Ulen_full,Pmap_full],)
"""

#注：只有LAT的93GHz和145GHz的noise map包含QUP noise(进行delensing), SAT只有QU noise(用于cross和fit).

"""
#FG patch correction
nus= {27: 27*1e9, 39 : 39*1e9, 93 : 93*1e9, 145 : 145*1e9, 225 : 225*1e9, 280 : 280*1e9,}
for idx in np.array([14,]):
    fg_map,_ = gaussian_fg(2048, nus, idx + 24325, same=False, method=1)
    hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/27/IQU_FG_{idx}.fits",fg_map[27],overwrite=True)
    hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/39/IQU_FG_{idx}.fits",fg_map[39],overwrite=True)
    hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/93/IQU_FG_{idx}.fits",fg_map[93],overwrite=True)
    hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/145/IQU_FG_{idx}.fits",fg_map[145],overwrite=True)
    hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/225/IQU_FG_{idx}.fits",fg_map[225],overwrite=True)
    hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/ALL_MAP_2048/FG/280/IQU_FG_{idx}.fits",fg_map[280],overwrite=True)
    """