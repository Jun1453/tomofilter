# This script applies the tomographic filters (R6000.out) to the numpy files for resampled models (xxx.mapped.npy)
# and then output filtered models (xxx.dag.npy).  Note that filters are n-by-n matrices where n is the length of
# resampled models. The .out files of HMSL-P06 and HMSL-S06 are ready on Zenodo (doi:10.5281/zenodo.6474359).

import numpy as np
import pandas as pd

class ResolutionOperator():
    def __init__(self, filename):
        data = np.fromfile(filename, dtype='i4')
        (M, N) = (data[0], data[1])
        del data
        self.matrix = np.fromfile(filename, dtype=np.float64)[1:].reshape(M, N, order='F')

class Model():
    def __init__(self, name, suffix=""):
        np.random.seed(int(id(self))%2**16)
        if "HMSL" in name:
            self.values = self._getanomalies(f"/Users/jun/mantletomo/legacy/ref{name[5]}Tomo.csv")
        else:
            self.values = np.load(f"{name}.mapped{suffix}.npy")
        self.filteredby = []
        self.modelname = name
    def __str__(self):
        return self.modelname
    def _getanomalies(self, filepath):
        table = pd.read_csv(filepath)
        return table['anomaly'].values
    def filter(self, filter: ResolutionOperator):
        self.values = filter.matrix @ self.values[0:filter.matrix.shape[0]]
        self.filteredby.append(filter)
        return self
        
Rp = ResolutionOperator('./Res_P/R6000.out') # revised R for VP
Rs = ResolutionOperator('./Res_S/R6000.out') # revised R for VS

##
# Filtering synthetic models
##

S_models = list(map(Model, ['T1_dlnVs', 'T1-pPv_dlnVs', 'TC1_dlnVs', 'TC1-pPv_dlnVs', 'TC4_dlnVs', 'TC4-pPv_dlnVs']))
P_models = list(map(Model, ['T1_dlnVp', 'T1-pPv_dlnVp', 'TC1_dlnVp', 'TC1-pPv_dlnVp', 'TC4_dlnVp', 'TC4-pPv_dlnVp']))

for model in S_models: 
    np.save(f"{model}.dag", model.filter(Rs).values)

for model in P_models: 
    np.save(f"{model}.dag", model.filter(Rp).values)

##
# Filtering HMSL models
##

np.save(f"StoP.dag", Model("HMSL-S06").filter(Rp).values)
np.save(f"PtoS.dag", Model("HMSL-P06").filter(Rs).values)
