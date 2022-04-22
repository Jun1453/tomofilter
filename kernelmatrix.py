def getSurfaceKernels(filenames):
    import numpy as np
    from scipy.io import FortranFile

    M = 2578 * len(filenames)
    N = 46404
    K = np.ndarray(shape=(M,N))
    residuals = []
    filenum = 0


    for filename in filenames:
        f = FortranFile(filename, 'r')
        count = 0
        while(1):
            if count % 2578 == 0 and count > 0:
                filenum += 1
                break
            else:
                mh = f.read_record(dtype='f4, f4, i4')
            if not mh: break
            (rem, res, nrow) = mh[0]
            
            rin = f.read_record(dtype=np.float32)
            kmat = f.read_record(dtype=np.int32)

            residuals.append(rem/(2*res))

            row = np.zeros(N)
            for j in range(nrow):
                row[kmat[j]-1] = rin[j] / (2*res)
            K[count + 2578*filenum] = row

            count += 1

    residuals = np.reshape(residuals,(len(residuals),1))
    return(residuals,K)

def getBodyKernel(fileprefix):
    import numpy as np
    import pandas as pd
    from scipy.io import FortranFile
    from math import sqrt

    table = pd.read_csv('%s.table' % (fileprefix), delim_whitespace=True)
    y = table['residual']
    error = table['error']
    cor1 = table['ellipcor']
    cor2 = table['crustcor']
    
    M = table.shape[0]
    N = 46404
    K = np.ndarray(shape=(M,N))
    f = FortranFile('%s.swp.mat' % (fileprefix), 'r')
    
    for i in range(M):
        mh = f.read_record(dtype='i4, a4, i4, i4')
        (nrow, sta, yr, dat) = mh[0] #pylint: disable=unused-argument
        rin = f.read_record(dtype=np.float32)
        kmat = f.read_record(dtype=np.int32)

        row = np.zeros(N)
        for j in range(nrow):
            row[kmat[j]-1] = rin[j] / error.values[i]
        K[i][:] = row

    d = (y.values[0:M]+cor1.values[0:M]+cor2.values[0:M])/error.values[0:M]
    d = np.reshape(d,(len(d),1))
    return(d,K)