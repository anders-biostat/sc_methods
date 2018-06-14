import numpy as np
cimport numpy as npc
#from cysignals.signals cimport sig_check


cdef npc.float64_t square( npc.float64_t x ):
   return x*x

cimport cython
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.

def dist( a, int nnz ):

    # sparse matrix format consists of 3 vectors
    indptr_np = np.zeros(a.shape[0] + 1, dtype=np.intc)
    indices_np = np.zeros(nnz, dtype=np.intc)
    data_np = np.zeros(nnz, dtype=np.float64)

    cdef int[:] a_indptr = a.indptr
    cdef int[:] a_indices = a.indices
    cdef npc.float64_t[:] a_data = a.data

    cdef int[:] indptr = indptr_np
    cdef int[:] indices = indices_np
    cdef npc.float64_t[:] data = data_np
 
    cdef int ncells = a.shape[0]

    cdef int indptr_i = 1
    cdef int data_i = 0

    cdef int cellA, cellB, indptrA, indptrB
    cdef npc.float64_t s, sim

    for cellA in range(ncells):
        for cellB in range(ncells):
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
                else:  # a.indices[indptrA] == a.indices[indptrB]
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
        #sig_check()
    return(data_np, indices_np, indptr_np)
