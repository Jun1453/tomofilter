# This script applies the tomographic filters (R6000.out) to the numpy files for resampled models (xxx.mapped.npy)
# and then output filtered models (xxx.dag.npy).  Note that filters are n-by-n matrices where n is the length of
# resampled models. The .out files of HMSL-P06 and HMSL-S06 are ready on Zenodo (doi:10.5281/zenodo.6474359).

import numpy as np
class ResolutionOperator():
    def __init__(self, filename):
        data = np.fromfile(filename, dtype='i4')
        (M, N) = (data[0], data[1])
        del data
        self.matrix = np.fromfile(filename, dtype=np.float64)[1:].reshape(M, N, order='F')

    def filter(self, model):
        mt = np.load(f"{model}.npy")
        return self.matrix @ mt[0:self.matrix.shape[0]]
        
Rp = ResolutionOperator('./Res_P/R6000.out') # revised R for VP
Rs = ResolutionOperator('./Res_S/R6000.out') # revised R for VS

# Filtering Synthetic models

P_models = ['T1_dlnVp', 'T1-pPv_dlnVp', 'TC1_dlnVp', 'TC1-pPv_dlnVp', 'TC4_dlnVp', 'TC4-pPv_dlnVp']
S_models = ['T1_dlnVs', 'T1-pPv_dlnVs', 'TC1_dlnVs', 'TC1-pPv_dlnVs', 'TC4_dlnVs', 'TC4-pPv_dlnVs']

for model in P_models: 
    np.save(f"{model}.dag", Rp.filter(f"{model}.mapped"))

for model in S_models: 
    np.save(f"{model}.dag", Rs.filter(f"{model}.mapped"))

# Cross-filtering
Stomo_model = "refStomo"
np.save(f"RpVs.dag", Rp.filter(Stomo_model))
