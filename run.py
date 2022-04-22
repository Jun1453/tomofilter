import subprocess
import h5sparse as h5py
from inversion import create_params, inversion
from matrix import createInputFile, read_matrices

filename = './database.hdf'

# implement svd  to create sensitivity matrix and placeholder singular matrix in hdf file
modelname = 'Sonly'
modelnumber = 100
params = create_params(f=modelname, k=modelnumber, c=True, z=True, smooth=False, over=True, l=2)
inversion(params, output_dir='./')
csr, d, c = read_matrices(filename, modelname, modelnumber)

## generate matrix G and d
createInputFile(csr.A, './kernel_input')
subprocess.run("./arrayconvert && \\rm -f kernel_input && mv double_input double_input.Gs", shell=True)
createInputFile(d, './kernel_input')
subprocess.run("./arrayconvert && \\rm -f kernel_input && mv double_input double_input.ds", shell=True)


# implement the same svd for P wave
modelname = 'P'
modelnumber = 100
params = create_params(f=modelname, k=modelnumber, c=True, z=True, smooth=False, over=True, l=2)
inversion(params, output_dir='./')
csr, d, c = read_matrices(filename, modelname, modelnumber)

## generate matrix G and d
createInputFile(csr.A, './kernel_input')
subprocess.run("./arrayconvert && \\rm -f kernel_input && mv double_input double_input.Gp", shell=True)
createInputFile(d, './kernel_input')
subprocess.run("./arrayconvert && \\rm -f kernel_input && mv double_input double_input.dp", shell=True)
