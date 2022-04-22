import config
import numpy as np
import pandas as pd
from math import sqrt
from scipy.io import FortranFile

h = list(map(lambda r1, r2: sqrt(r2-r1), config.RADIUS1, config.RADIUS2))

def body_matrix(fileprefix):
    table = pd.read_csv(f"{config.MATPATH}/{fileprefix}.table", delim_whitespace=True)
    M = table.shape[0]
    y = table['residual']
    error = table['error']
    cor1 = table['ellipcor']
    cor2 = table['crustcor']
    
    K = np.ndarray(shape=(M,config.N))
    f = FortranFile(f"{config.MATPATH}/{fileprefix}.swp.mat", 'r')
    
    for i in range(M):
        mh = f.read_record(dtype='i4, a4, i4, i4')
        (nrow, sta, yr, dat) = mh[0]
        rin = f.read_record(dtype=np.float32)
        kmat = f.read_record(dtype=np.int32)

        row = np.zeros(config.N)
        for j in range(nrow):
            row[kmat[j]-1] = rin[j] / error.values[i]
            # if is_thickness_to_weight: row[kmat[j]-1] /= h[(kmat[j]-1)//2578]
        K[i][:] = row
        # if i == 1:
            # print(csr_matrix(row).nonzero())

    residuals = (y.values[0:M]+cor1.values[0:M]+cor2.values[0:M])/error.values[0:M]
    residuals = np.reshape(residuals,(len(residuals),1))
    return residuals, K

def surface_matrix(filenames):
    K = np.ndarray(shape=(2578*len(filenames), config.N))
    residuals = []

    for i in range(len(filenames)):
        f = FortranFile(f"{config.MATPATH}/{filenames[i]}", 'r')
        count = 0
        while(1):
            if count == 2578: break
            mh = f.read_record(dtype='f4, f4, i4')
            if not mh: break
            (rem, res, nrow) = mh[0]
            
            rin = f.read_record(dtype=np.float32)
            kmat = f.read_record(dtype=np.int32)

            residuals.append(rem/(2*res))

            row = np.zeros(config.N)
            for j in range(len(kmat)):
            # for j in range(nrow):
                row[kmat[j]-1] = rin[j] / (2*res)
                # if is_thickness_to_weight: row[kmat[j]-1] /= h[(kmat[j]-1)//2578] 
            K[count + 2578*i] = row

            count += 1

    residuals = np.reshape(residuals,(len(residuals),1))
    return residuals, K

def revert_h_weighting(m):
    newm = [ m[j]/h[j//2578] for j in range(len(m))]
    return np.array(newm)

def revert_c_weighting(m, c):
    return m[:,0]/c if len(m.shape)==2 else m[:]/c

def createInputFile(A, filename):

    f = FortranFile(filename, 'w')
    if len(A.shape) == 2:
        (M, N) = A.shape
    else:
        (M, N) = (len(A), 1)
    f.write_record(int(M))
    f.write_record(int(N))
    # f.write_record(np.array([M,N], dtype=int))
    # for val in A.flatten(order='F'):
    #     f.write_record(float(val))
    for i in range(N):
        this_column = A[:,i] if N>1 else A[:]
        this_column = np.array(this_column, dtype=np.float64)
        f.write_record(this_column)

def readOutputFile(filename):
    data = np.fromfile(filename, dtype='i4')
    (M, N) = (data[0], data[1])
    del data
    data = np.fromfile(filename, dtype=np.float64)[1:].reshape(M, N, order='F')
    return data