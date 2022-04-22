import os
import numpy as np
import h5sparse as h5py
import subprocess
from scipy.io import FortranFile

def createInputFile(A, filename):
    f = FortranFile(filename, 'w')
    if len(A.shape) == 2:
        (M, N) = A.shape
    else:
        (M, N) = (len(A), 1)
    f.write_record(int(M))
    f.write_record(int(N))
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
    
def open_hdf5(filename):
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    try:
        db = h5py.File(filename, mode='a')
    except:
        db = h5py.File(filename, mode='w')
    return db

def read_matrices(filename, modelname, modelnumber):
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    file = h5py.File(filename, mode='r')
    csr = file[f"/{modelname}/{modelnumber}/kernel"][()]
    d = file[f"/{modelname}/{modelnumber}/residual"][...]
    c = file[f"/{modelname}/{modelnumber}/weight"][...]
    file.close()
    return csr, d, c

# run as main
# if __name__ == '__main__':
#     ## load matrix data
#     filename = './database.hdf'
#     f = h5py.File(filename, mode='r')
#     csr = f["/Sonly/100/kernel"][()]
#     d = f["/Sonly/100/residual"][...]
#     c = f["/Sonly/100/weight"][...]
#     f.close()

#     ## prepare the matrices to solve the equation Gm = d
#     ## generate matrix G
#     createInputFile(csr.A, './kernel_input')
#     subprocess.run("./arrayconvert && mv double_input double_input.G", shell=True)

#     ## generate matrix d
#     createInputFile(d, './kernel_input')
#     subprocess.run("./arrayconvert && mv double_input double_input.d", shell=True)

