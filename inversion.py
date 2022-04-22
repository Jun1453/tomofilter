import argparse
import h5sparse as h5py
from h5py import SoftLink
import os
import numpy as np
import blockmodel
from math import sqrt
from read_matrix import body_matrix, surface_matrix, revert_h_weighting, revert_c_weighting
from scipy.sparse import csr_matrix, diags, vstack
from scipy.sparse.linalg import svds
import pandas as pd
from time import time
from datetime import date

# kernel object
class Kernel():
    ## initialized with filegroup name & k value
    def __init__(self, params):
        self.filegroup = params.f
        self.num_of_vectors = params.k
        self.epsilon = params.epsilon
        self.regularization = params.regular
        self.is_coverage_to_weight = params.c
        self.is_thickness_to_weight = params.z
        self.is_smoother_to_apply = params.smooth
        self.is_loaded = False
        self.smoother = False

    ## load kernel results
    def load_kernel(self, db):
        if f"/{self.filegroup}" in db:
            for item in db[f"/{self.filegroup}"]:
                if f"/{self.filegroup}/{item}/kernel" in db:
                    attrs = db[f"/{self.filegroup}/{item}/kernel"].attrs
                    if attrs['thickness_weighted'] == self.is_thickness_to_weight:
                        if (not self.is_coverage_to_weight and attrs['coverage_weighted'] == -1) or (self.is_coverage_to_weight and attrs['coverage_weighted'] == self.regularization):
                            if not 'smoother_added' in attrs: attrs['smoother_added'] = False
                            if attrs['smoother_added'] == self.is_smoother_to_apply:
                                self.csr = db[f"/{self.filegroup}/{item}/kernel"][()]
                                self.d = db[f"/{self.filegroup}/{item}/residual"]
                                self.smoother = attrs['smoother_added']
                                if self.is_coverage_to_weight: self.c = db[f"/{self.filegroup}/{item}/weight"]
                                if not item == self.num_of_vectors:
                                    if not f"/{self.filegroup}/{self.num_of_vectors}/kernel" in db: db[f"/{self.filegroup}/{self.num_of_vectors}/kernel"] = SoftLink(f"/{self.filegroup}/{item}/kernel")
                                    if not f"/{self.filegroup}/{self.num_of_vectors}/residual" in db: db[f"/{self.filegroup}/{self.num_of_vectors}/residual"] = SoftLink(f"/{self.filegroup}/{item}/residual")
                                    if not f"/{self.filegroup}/{self.num_of_vectors}/weight" in db and self.is_coverage_to_weight: db[f"/{self.filegroup}/{self.num_of_vectors}/weight"] = SoftLink(f"/{self.filegroup}/{item}/weight")
                                self.is_loaded = True
                                break
    
    ## save kernel results
    def save_kernel(self, db):
        if not f"/{self.filegroup}" in db: db.create_group(f"/{self.filegroup}")
        if not f"/{self.filegroup}/{self.num_of_vectors}" in db: db.create_group(f"/{self.filegroup}/{self.num_of_vectors}")

        kernel_data = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/kernel", self.csr)
        kernel_data.attrs['created'] = date.today().strftime('%Y-%m-%d')
        kernel_data.attrs['error_weighted'] = True
        kernel_data.attrs['smoother_added'] = self.smoother
        kernel_data.attrs['thickness_weighted'] = self.is_thickness_to_weight
        kernel_data.attrs['coverage_weighted'] = self.regularization if self.is_coverage_to_weight else -1

        if self.is_coverage_to_weight:
            weight_data = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/weight", self.c)
            weight_data.attrs['created'] = date.today().strftime('%Y-%m-%d')
            weight_data.attrs['smoother_added'] = self.smoother

        residual_data = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/residual", self.d)
        residual_data.attrs['created'] = date.today().strftime('%Y-%m-%d')
        residual_data.attrs['error_weighted'] = True
        residual_data.attrs['smoother_added'] = self.smoother

    ## Load existing data matrix
    def readmatrix(self):
        if self.filegroup == 'Sfast':
            filenames = ['ScS-S.4.05', 'SS-S.4']
            is_surface_to_add = True
        elif self.filegroup == 'Surf':
            filenames = []
            is_surface_to_add = True
        elif self.filegroup == 'Scomb':
            filenames = ['Scomb.4']
            is_surface_to_add = False
        elif self.filegroup == 'ScombSurf':
            filenames = ['Scomb.4']
            is_surface_to_add = True
        elif self.filegroup == 'Sonly':
            filenames = ['Scomb.4', 'SS.95-05.4', 'ScS-S.4.05', 'SS-S.4']
            is_surface_to_add = False
        elif self.filegroup == 'S':
            filenames = ['Scomb.4', 'SS.95-05.4', 'ScS-S.4.05', 'SS-S.4']
            is_surface_to_add = True
        elif self.filegroup == 'P':
            filenames = ['Pcomb.4', 'PP.95-05.4',  'PP-P.4']
            is_surface_to_add = False
        else:
            raise ValueError('unknown file combination %s' % self.filegroup)

        data = ()
        kernels = ()
        for filename in filenames:
            (d, K) = body_matrix(filename)
            data = data + (d,)
            kernels = kernels + (K,)
        if is_surface_to_add:
            fns = ['mat4.nd.r'+str(i)+'.swp' for i in range(4,15+1)]
            fns += ['mat4.nd.l'+str(i)+'.swp' for i in range(4,15+1)]
            (d, K) = surface_matrix(fns)
            data = data + (d,)
            kernels = kernels + (K,)

        return np.vstack(data), np.vstack(kernels)

    ## collect and weight kernel matrices
    def get_kernel(self):
        from config import H
        if not self.is_thickness_to_weight:
            H = np.ones(len(H))
        (self.d, K) = self.readmatrix()

        # Apply coverage weighting column by column (save memory)
        c = np.zeros(K.shape[1])
        for i in range(len(c)):
            c[i] = H[i//2578]
        if self.is_coverage_to_weight:
            for i in range(len(c)):
                if self.regularization == 0:
                    c[i] *= max(1, np.sum(K[:,i]!=0))
                elif self.regularization == 1:
                    c[i] *= max(1, np.sum(np.abs(K[:,i]), axis=0))
                elif self.regularization == 2:
                    c[i] *= max(1, np.sqrt(np.sum(K[:,i]**self.regularization, axis=0)))
                elif self.regularization > 2:
                    c[i] *= max(1, np.sqrt(np.sum(np.abs(K[:,i])**self.regularization, axis=0)))
                else:
                    raise ValueError('incorrect regularization norm: L-%s' % self.regularization)
                K[:,i] = K[:,i]/c[i]
        self.c = c

        # Compress K into scipy sparse matrix (csr format)
        self.csr = csr_matrix(K,dtype=np.float32)
        self.is_loaded = True

    ## Apply first differnce smoother
    def first_diffence(self, db, radial_lambda=5., lateral_lambda=5.):
        if self.smoother: return
        m = blockmodel.Model()
        dwt = self.c[...] ** 2

        # Insert extra rows to sensitivity matrix
        data = []
        row = []
        col = []

        # Insert extra elements to traveltime residual
        (mshell, nshell) = (2578, 18)

        # Radical smoother
        for i in range(mshell):
            for j in range(nshell-1):
                lower_block_num = i + j * mshell
                upper_block_num = lower_block_num + mshell

                # insert flam for lower block and -flam for upper to matrix
                newrow = row[-1] + 1 if len(row) > 0 else 0
                row.append(newrow)
                col.append(lower_block_num)
                data.append(radial_lambda / self.c[lower_block_num])
                dwt[lower_block_num] += (radial_lambda / self.c[lower_block_num]) ** 2

                row.append(newrow)
                col.append(upper_block_num)
                data.append(-radial_lambda / self.c[upper_block_num])
                dwt[upper_block_num] += (radial_lambda / self.c[upper_block_num]) ** 2

                # insert 0 residual to self.d
                self.d = np.append(self.d, 0)
                
        # Lateral smoother
        for this_block in range(nshell*mshell):
            # get surrounding blocks
            west_neigbor = m[this_block].neighbor('W').id
            east_neigbor = m[this_block].neighbor('E').id
            south_neighbor = m[this_block].neighbor('S').id
            north_neighbor = m[this_block].neighbor('N').id

            # insert 4*glam for center block and -glam for ambient to matrix
            newrow = row[-1] + 1 if len(row) > 0 else 0
            row.append(newrow)
            col.append(this_block)
            data.append(4 * lateral_lambda / self.c[this_block])
            dwt[this_block] += (4 * lateral_lambda / self.c[this_block]) ** 2

            row.append(newrow)
            col.append(north_neighbor)
            data.append(-lateral_lambda / self.c[north_neighbor])
            dwt[north_neighbor] += (lateral_lambda / self.c[north_neighbor]) ** 2

            row.append(newrow)
            col.append(south_neighbor)
            data.append(-lateral_lambda / self.c[south_neighbor])
            dwt[south_neighbor] += (lateral_lambda / self.c[south_neighbor]) ** 2

            row.append(newrow)
            col.append(east_neigbor)
            data.append(-lateral_lambda / self.c[east_neigbor])
            dwt[east_neigbor] += (lateral_lambda / self.c[east_neigbor]) ** 2

            row.append(newrow)
            col.append(west_neigbor)
            data.append(-lateral_lambda / self.c[west_neigbor])
            dwt[west_neigbor] += (lateral_lambda / self.c[west_neigbor]) ** 2

            # insert 0 residual to self.d
            self.d = np.append(self.d, 0)

        # compile smoother matrix and stack with csr
        print(f"inserting {len(row)} elements smoother in {row[-1]+1} rows...")
        smoother = csr_matrix((data, (row, col)), shape=(row[-1]+1, mshell*nshell))
        self.csr = vstack([self.csr, smoother])
        self.c = np.sqrt(dwt)

        # toggle on smoother status
        self.smoother = True

    ## decompose and reconstruct
    def svd(self, db):
        t = time()
        U, W, VT = svds(self.csr, k=self.num_of_vectors) if self.num_of_vectors < 46404 else svds(self.csr, k=46403)
        print('SVD with %s vectors finished in %.2f s' % (str(self.num_of_vectors) if self.num_of_vectors < 46404 else 'All', time()-t))

        ### Reshape and save SVD results with timing
        t = time()
        thres = self.epsilon * max(W)
        z_array = []
        for i in range(W.shape[0]):
            if W[i] < thres:
                break
            z_array.append(1/W[i])
        Z = diags(z_array)
        rank = len(z_array)

        V = np.transpose(VT[0:rank, :])
        UT = np.transpose(U[:, 0:rank])

        ### write to database
        svd_v_data = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/svd_v", V)
        svd_z_data = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/svd_z", z_array)
        svd_ut_data = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/svd_ut", UT)
        svd_z_data.attrs['epsilon'] = self.epsilon

        t = time()
        m = V @ (Z @ (UT @ self.d))
        if self.is_thickness_to_weight: m = revert_h_weighting(m)
        if self.is_coverage_to_weight: m = revert_c_weighting(m, self.c)
        s = pd.Series(m.reshape(len(m)))
        model = overwrite_data(db, f"/{self.filegroup}/{self.num_of_vectors}/model", m)
        print(s.describe())
        print('inversion model calculated from SVD in %.2f s' % (time()-t))


# parse input parameters
def create_parser():
    parser = argparse.ArgumentParser(prog="Tomography Inversion Tools")
    parser.add_argument("f", help="filegroup")
    parser.add_argument("-k", help="truncated size", type=int, default=6000)
    parser.add_argument("--epsilon", "-e", help="neglecting threshold", type=float, default=1e-4)
    parser.add_argument("--regular", "-l", help="power magnitude for regularization", type=int, default=2)
    parser.add_argument("-c", help="coverage weighting", action='store_true')
    parser.add_argument("-z", help="thickness weighting", action='store_true')
    parser.add_argument("--smooth", help="apply smoother", action='store_true')
    parser.add_argument("--overwrite", "-o", help="overwrite kernel", action='store_true')
    parser.add_argument("--dispose", "-d", help="not to save kernel", action='store_true')
    return parser

def create_params(f, **kwargs):
    parser = create_parser()
    args = []
    args.append(f)
    for key in kwargs.keys():
        if not kwargs[key] is False:
            args.append(f"--{key}" if len(key) > 1 else f"-{key}") # note this is only compatible with Python >= 3.6
            if not type(kwargs[key]) is bool:args.append(str(kwargs[key]))
    return parser.parse_args(args)

# main program
def inversion(params, output_dir='.'):
    ## prepare kernel data
    db = open_hdf5(f"{output_dir}/database.hdf")
    data = Kernel(params)

    if not params.overwrite:
        print("try to load sensitivity matrix if possible...")
        data.load_kernel(db)

    if not data.is_loaded:
        if not params.overwrite: print("no sensitivity matrix availible")
        print("generating sensitivity matrix from data...")
        data.get_kernel()
        ## add smoother to kernel matrix
        if params.smooth: data.first_diffence(db, radial_lambda=5, lateral_lambda=5)
        if not params.dispose:
            data.save_kernel(db)
            print("sensitivity matrix saved.")
    else:
        print("sensitivity matrix loaded.")

    print(f"sensitivity matrix shape: {data.csr.shape}")
    ## implement svd
    data.svd(db)
    db.close()
    
# open database
def open_hdf5(filename):
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    try:
        db = h5py.File(filename, mode='a')
    except:
        db = h5py.File(filename, mode='w')
    return db

# handle database write
def overwrite_data(db, path, newdata):
    if path in db: del db[path]
    data = db.create_dataset(path, data=newdata)
    return data


# run as main
if __name__ == '__main__':
    parser = create_parser()
    params = parser.parse_args()
    inversion(params)
