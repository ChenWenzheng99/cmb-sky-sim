from mpi4py import MPI
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

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
sys.path.append('/sharefs/alicpt/users/chenwz/Testarea/FORECAST/Foreground/foreground_run')

from lensingb_mine import *
from noise_making import *
from utils_mine import *
from phi_noise_making import *
from generate_map import *
from lenspyx.utils import camb_clfile
from lenspyx.utils_hp import synalm, almxfl, alm2cl

camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))    
sys.path.insert(0,camb_path)

import camb
from camb import model, initialpower ,correlations

from scipy.special import factorial

import pymaster as nm


# 获取传递给脚本的参数
args = sys.argv[1:]

# 打印当前作业的编号
# 计算当前作业的文件编号列表
file_indices = comm_rank + int(args[0])

# 打印当前作业的文件编号列表
print(f"Job {comm_rank} will process files: {file_indices}")


def kappa_to_phi(plm, lmin, lmax, wiener=None):
    q2k = lambda l: l*(l + 1) / 2 # potential -> convergence
    q2d = lambda l: (l*(l + 1)) ** 0.5 # potential -> deflection

    dlmax = hp.Alm.getlmax(np.shape(plm)[0])
    cut = np.where((np.arange(dlmax+1 ) >= lmin) * (np.arange(dlmax+1) <= lmax), 1, 0) # band limit
    if wiener is None:
        wiener = np.ones(dlmax + 1)
    return hp.almxfl(plm, wiener * cut * cli(q2k(np.arange(dlmax + 1))))

from plancklens import utils as utils_pk
cl_len = utils_pk.camb_clfile('/sharefs/alicpt/users/chenwz/download/plancklens/plancklens/data/cls/FFP10_wdipole_lensedCls.dat')
cl_unl = utils_pk.camb_clfile('/sharefs/alicpt/users/chenwz/download/plancklens/plancklens/data/cls/FFP10_wdipole_lenspotentialCls.dat')
ell = np.arange(len(cl_unl['tt']))

lmax_unl = mmax_unl = 3*2048 - 1 - 144
dlmax = 3*2048 - 1

alm_unl = hp.synalm([cl_unl['tt'], cl_unl['te'], cl_unl['ee'], cl_unl['bb']], lmax=lmax_unl, mmax=mmax_unl)


kappa_2048 = hp.read_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048/{file_indices:05d}/kappa_{file_indices:05d}.fits", verbose=False)
plm = kappa_to_phi(hp.map2alm(kappa_2048, lmax=3*2048-1), 0, 3*2048-1, wiener=None)
dlm = hp.almxfl(plm, np.sqrt(np.arange(dlmax + 1, dtype=float) * np.arange(1, dlmax + 2)), None, False)  

geom_info = ('healpix', {'nside':2048})
Tlen, Qlen, Ulen = lenspyx.alm2lenmap(alm_unl, dlm, geometry=geom_info, verbose=1, epsilon=1e-6)

hp.write_map(f"/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048/{file_indices:05d}/lensedcmb_new_{file_indices:05d}.fits", [Tlen, Qlen, Ulen], overwrite=True)