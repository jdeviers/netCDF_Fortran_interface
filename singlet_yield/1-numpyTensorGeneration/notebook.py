import math
import functools
import numpy as np
import scipy as sci
from scipy import linalg
import matplotlib.pyplot as plt
import qutip as qt
from numba import jit, prange, float64, complex128, int64

#from py3_f90_evalYield import evalyield as f90_mod

opstr2fun = {'x': lambda dim: qt.spin_Jx((dim-1)/2),
             'y': lambda dim: qt.spin_Jy((dim-1)/2),
             'z': lambda dim: qt.spin_Jz((dim-1)/2),
             'p': lambda dim: qt.spin_Jp((dim-1)/2),
             'm': lambda dim: qt.spin_Jm((dim-1)/2),
             'i': qt.identity}

def mkSpinOp(dims, specs):
    ops = [qt.identity(d) for d in dims]
    for ind, opstr in specs:
        ops[ind] = ops[ind] * opstr2fun[opstr](dims[ind])
    return qt.tensor(ops)
def mkH1(dims, ind, parvec):
    axes = ['x', 'y', 'z']
    return functools.reduce(lambda a, b: a + b, 
               [v * mkSpinOp(dims, [(ind,ax)]) for v, ax in zip(parvec, axes) if v!=0])
def mkH12(dims, ind1, ind2, parmat):
    axes = ['x', 'y', 'z']
    ops = []
    for i in range(3):
        for j in range(3):
            if parmat[i,j] != 0:
                ops.append(parmat[i,j] * mkSpinOp(dims, [(ind1,axes[i]), (ind2,axes[j])]))
    return functools.reduce(lambda a, b: a + b, ops)

# FAD*- hfcs in MHz
hfcsF = {
    'N5' : [[-2.84803, 0.0739994, -1.75741],
            [0.0739994, -2.5667, 0.326813],
            [-1.75741, 0.326813, 53.686]],
    'N10': [[-0.0979402, 0.00195169, 1.80443],
            [0.00195169, -0.513124, -0.508695],
            [1.80443, -0.508695, 19.109]],
    'H7' : [[-3.7789, 0., 0.],
            [0., -3.7789, 0.],
            [0., 0., -3.7789]],
    'H8' : [[11.872, 0., 0.],
            [0., 11.872, 0.],
            [0., 0., 11.872]],
    'Hb1': [[8.48607, -0.878249, -1.22263],
            [-0.878249, 5.56589, 0.299019],
            [-1.22263, 0.299019, 5.3476]],
    'Hb2': [[5.30097, 1.02387, -1.09139],
            [1.02387, 2.32064, -0.269562],
            [-1.09139, -0.269562, 1.93875]],
    'H9' : [[2.3771, 0.860593, -0.0221318],
            [0.860593, 4.16473, 0.177746],
            [-0.0221318, 0.177746, 0.324792]],
    'H6' : [[-5.31431, -1.03422, 0.256356],
            [-1.03422, -13.3929, -0.200838],
            [0.256356, -0.200838, -11.6316]],
    'Hc' : [[1.27891, -0.0243844, 0.144687],
            [-0.0243844, -0.595294, -0.014526],
            [0.144687, -0.014526, -0.568057]]
    }
# Remarks: odd number of H -> <m|Sx,y,z|m> = 0 -> diagonal contrib = 1/4
#          repeat H7 and H8 if desired
# nucLabelsF = ['N5', 'N10', 'H6', 'H8', 'Hb1', 'Hb2', 'H9', 'H7', 'Hc']
nucLabelsF = ['N5', 'N10', 'H6']
As1 = [np.array(hfcsF[nuc])*2*math.pi for nuc in nucLabelsF]
dims1 = [2, *[3 if nuc[0]=='N' else 2 for nuc in nucLabelsF]]

print(dims1)
print(len(As1))
for i in range(len(As1)):
    print(i)
    print( mkH12(dims1, 0, i+1, As1[i]) )
    
print(mkSpinOp(dims1, [(0,'x')]))
    
H1 = sum(mkH12(dims1, 0, i+1, As1[i]) for i in range(len(As1)))
lambda1, V = np.linalg.eigh(H1.full())
Sxyz1 = np.array([V.conj().T @ (mkSpinOp(dims1, [(0,ax)]).full() @ V) for ax in ['x', 'y', 'z']])

# W*+
hfcsW = {
    'N1' : [[-1.94218, -0.0549954, -0.21326],
            [-0.0549954, -2.29723, -0.441875],
            [-0.21326, -0.441875, 19.156]],
    'H1' : [[-2.14056, 6.31534, 0.17339],
            [6.31534, -18.9038, -0.0420204],
            [0.17339, -0.0420204, -14.746]],
    'H2' : [[-21.1751, 4.41952, 0.163566],
            [4.41952, -4.32747, 0.110325],
            [0.163566, 0.110325, -15.993]],
    'Hb1': [[8.39562, -2.71765, -0.582406],
            [-2.71765, 8.11649, 0.637991],
            [-0.582406, 0.637991, 4.94689]],
    'Hb2': [[27.4878, 0.814461, -1.79339],
            [0.814461, 24.0013, -0.351985],
            [-1.79339, -0.351985, 22.8546]],
    'H4' : [[-6.44205, 0.787534, 0.270961],
            [0.787534, -23.0383, -0.0378865],
            [0.270961, -0.0378865, -17.0576]],
    'H5' : [[2.50597, -1.01627, 0.0195079],
            [-1.01627, 4.51584, -0.0838277],
            [0.0195079, -0.0838277, 0.654643]],
    'H6' : [[-14.5993, -4.81629, 0.0840456],
            [-4.81629, -5.12283, 0.0246775],
            [0.0840456, 0.0246775, -11.2632]],
    'H7' : [[-1.05186, 1.16764, 0.0917152],
            [1.16764, -8.50485, 0.0678854],
            [0.0917152, 0.0678854, -7.36221]]
    }
# nucLabelsW= ['N1', 'H1', 'H2', 'Hb1', 'Hb2', 'H4', 'H6', 'H7', 'H5']
nucLabelsW= ['N1', 'H1']
As2 = [np.array(hfcsW[nuc])*2*math.pi for nuc in nucLabelsW]
dims2 = [2, *[3 if nuc[0]=='N' else 2 for nuc in nucLabelsW]]
H2 = sum(mkH12(dims2, 0, i+1, As2[i]) for i in range(len(As2)))
lambda2, V = np.linalg.eigh(H2.full())
Sxyz2 = np.array([V.conj().T @ (mkSpinOp(dims2, [(0,ax)]).full() @ V) for ax in ['x', 'y', 'z']])

Sxyz1T = np.transpose(Sxyz1, (1,2,0)).copy()
Sxyz2T = np.transpose(Sxyz2, (1,2,0)).copy()

Sxyz1T.flags


# ---------- Write Sxyz arrays to netCDF files ----------
import netCDF4 as nc

#try: ncfile.close()
#except: pass

ncfile = nc.Dataset('../new.nc', mode = 'w', format = 'NETCDF4')

ncfile.title = 'SxyzT 1 and 2'
ncfile.subtitle = 'A subtitle'

#print(np.amax(Sxyz1T)) # OK: same as read by f90 prog
#print(np.amin(Sxyz1T)) # OK.

#np.info(Sxyz1T) ; np.info(lambda1)
#np.info(Sxyz2T) ; np.info(lambda2)

# Setting dimensions for Sxyz1T
dim_S11 = ncfile.createDimension('S1_axis_1', Sxyz1T.shape[0])
dim_S12 = ncfile.createDimension('S1_axis_2', Sxyz1T.shape[1])

# Setting dimensions for Sxyz2T
dim_S21 = ncfile.createDimension('S2_axis_1', Sxyz2T.shape[0])
dim_S22 = ncfile.createDimension('S2_axis_2', Sxyz2T.shape[1])

# Setting dimensions for both tensors
dim_xyz = ncfile.createDimension('S_xyz', 3)           # 
dim_ReIm = ncfile.createDimension('Real_Imaginary', 2) # There is no complex type; saving as two REAL tensors

for dim in ncfile.dimensions.items():
	print(dim)

# Create the variables in nc to hold the tensors
nc_Sxyz1T = ncfile.createVariable('Sxyz1T', 'f8', ('S1_axis_1','S1_axis_2','S_xyz','Real_Imaginary'))
nc_Sxyz2T = ncfile.createVariable('Sxyz2T', 'f8', ('S2_axis_1','S2_axis_2','S_xyz','Real_Imaginary'))
nc_lambda1 = ncfile.createVariable('lambda1', 'f8', ('S1_axis_1',))
nc_lambda2 = ncfile.createVariable('lambda2', 'f8', ('S2_axis_1',))

# Fill the Sxyz1T var:
nc_Sxyz1T[:,:,:,0] = np.real(Sxyz1T) # With the Real part
nc_Sxyz1T[:,:,:,1] = np.imag(Sxyz1T) # and the Imaginary

# Same with Sxyz2T
nc_Sxyz2T[:,:,:,0] = np.real(Sxyz2T)
nc_Sxyz2T[:,:,:,1] = np.imag(Sxyz2T)

# And with the lambdas
nc_lambda1[:] = lambda1
nc_lambda2[:] = lambda2

print(ncfile)
ncfile.close()
# -------------------------------------------------------


@jit(nopython=True)
def evalYield_diag(k, Sxyz1, lambda1, Sxyz2, lambda2):
    # Sxyz = [Sx, Sy, Sz]
    # Sxyz1 = np.transpose(Sxyz1, (1,2,0))
    # Sxyz2 = np.transpose(Sxyz2, (1,2,0))
    d1 = Sxyz1.shape[0]
    d2 = Sxyz2.shape[0]
    z = d1 * d2 // 4
    v = 0.0
    for a1 in range(d1):
        sA = Sxyz1[a1,a1,:]
        for b1 in range(d2):
            sB = Sxyz2[b1,b1,:]
            #print(np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2)
            v += np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2
    v /= z
    return 1/4 + v

@jit(nopython=True,fastmath=True)
def evalYield_offdiag(k, Sxyz1, lambda1, Sxyz2, lambda2):
    # Sxyz = [Sx, Sy, Sz]
    # Sxyz1 = np.transpose(Sxyz1, (1,2,0))
    # Sxyz2 = np.transpose(Sxyz2, (1,2,0))
    d1 = Sxyz1.shape[0]
    d2 = Sxyz2.shape[0]
    z = d1 * d2 // 4
    v = 0.0
    k2 = k*k
    for a1 in range(d1):
        lambda1_a1 = lambda1[a1]
        for b1 in range(d2):
            a2 = a1
            b2 = b1
            sA = Sxyz1[a1,a2,:]
            dl1 = lambda1_a1 - lambda1[a2]
            while True:
                b2 += 1
                if b2 == d2:
                    b2 = 0
                    a2 += 1
                    if a2 == d1:
                        break
                    sA = Sxyz1[a1,a2,:]
                    dl1 = lambda1_a1 - lambda1[a2]
                sB = Sxyz2[b1,b2,:]
                dl2 = lambda2[b1] - lambda2[b2]
                v += np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2 / (k2 + (dl1 + dl2)**2)
    v *= k2/z*2
    return v

#evalYield_diag(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2)

#evalYield_diag(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2) + evalYield_offdiag(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2)

import random

@jit(nopython=True)
def evalYield_offDiag_random(k, Sxyz1, lambda1, Sxyz2, lambda2, nr_draws):
    # Sxyz = [Sx, Sy, Sz]
    # Sxyz1 = np.transpose(Sxyz1, (1,2,0))
    # Sxyz2 = np.transpose(Sxyz2, (1,2,0))
    d1 = Sxyz1.shape[0]
    d2 = Sxyz2.shape[0]
    z = d1 * d2 // 4
    v = 0.0
    k2 = k*k
    n = 0
    while n < nr_draws:
        a1 = np.random.randint(0, d1)
        a2 = np.random.randint(0, d1)
        b1 = np.random.randint(0, d2)
        b2 = np.random.randint(0, d2)
        if (a1 == a2) and (b1 == b2):
            continue
        dl1 = lambda1[a1] - lambda1[a2]
        dl2 = lambda2[b1] - lambda2[b2]
        sA = Sxyz1[a1,a2,:]
        sB = Sxyz2[b1,b2,:]
        v += np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2 / (k2 + (dl1 + dl2)**2)
        n += 1
    return (v / n) * k2/z

d1 = Sxyz1T.shape[0]
d2 = Sxyz2T.shape[0]
d = d1*d2
#(d**2-d)*evalYield_offDiag_random(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2, 1000000) + evalYield_diag(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2)

#(d**2.-d) # scaling could be included in function

#logn = np.arange(2.5,8,0.5)
#y0 = evalYield_offdiag(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2)
#print(y0)
#y = [(d**2.-d)*evalYield_offDiag_random(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2, n) for n in 10**logn]
#print(y)

#logn = np.arange(2.5,8,0.5)
#x0 = f90_mod.evalyield_offdiag2p(1.0, Sxyz1, lambda1, Sxyz2, lambda2)
#print(x0)
#x = [(d**2.-d)*f90_mod.evalyield_offdiag_random_omp(1.0, Sxyz1T, lambda1, Sxyz2T, lambda2, n) for n in 10**logn]
#print(x)

