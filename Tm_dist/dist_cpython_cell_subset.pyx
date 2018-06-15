cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np
import scipy as sp

ctypedef np.float32_t myfloat
ctypedef np.int32_t myint


cdef myfloat square(myfloat x):
    return x*x


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cdist( a, cell_subset = None ):
    '''
    a: a sparse matrix in csr format (float32)
    cell_subset: a list (or 1D numpy array) of cell indices;
        all cells will be used if no cell_subset is specified
    '''
    if cell_subset is None:
        cell_subset = range(a.shape[0])
    cdef myint[:] cell_sub = np.array(cell_subset, dtype=np.int32)

    cdef myint indptr_i, data_i, indptrA, indptrB, cellA, cellB
    cdef myfloat s, sim
    cdef myint n_cells = a.shape[0]
    cdef myint nnz = a.shape[0] * a.shape[0] / 5

    cdef myint[:] a_indptr = a.indptr
    cdef myint[:] a_indices = a.indices
    cdef myfloat[:] a_data = a.data
    
    # sparse matrix format consists of 3 vectors
    indptr_np = np.zeros(a.shape[0] + 1, dtype=np.int32)
    indices_np = np.zeros(nnz, dtype=np.int32)
    data_np = np.zeros(nnz, dtype=np.float32)
    cdef myint[:] indptr = indptr_np
    cdef myint[:] indices = indices_np
    cdef myfloat[:] data = data_np

    indptr_i = 1
    data_i = 0

    for cellA in range(n_cells):
        for cellB in cell_sub:
            indptrA = a_indptr[cellA]
            indptrB = a_indptr[cellB]
            s = 0
            while indptrA < a_indptr[cellA+1] and indptrB < a_indptr[cellB+1]:
                if a_indices[indptrA] < a_indices[indptrB]:
                    s += square( a_data[indptrA] )
                    indptrA += 1
                elif a_indices[indptrA] > a_indices[indptrB]:
                    s += square( a_data[indptrB] )
                    indptrB += 1
                else:  # a_indices[indptrA] == a_indices[indptrB]
                    s += square( a_data[indptrA] - a_data[indptrB] )
                    indptrA += 1
                    indptrB += 1
                if s > 0.4:
                    break
            else:
                sim = 1 - s/2
                # write output to sparse matrix (3 vectors):
                data[data_i] = sim
                indices[data_i] = cellB
                data_i += 1
                if data_i >= nnz:
                    raise IndexError( "More nonzeroes than specified as 'nnz'." )
        indptr[indptr_i] = data_i
        indptr_i += 1
    return sp.sparse.csr_matrix( (data_np, indices_np, indptr_np),
                                  shape = (n_cells, n_cells)).transpose()
