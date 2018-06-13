def dist(a):

    # sparse matrix format consists of 3 vectors
    indptr = np.zeros(a.shape[0] + 1, dtype="int32")
    indices = np.zeros(nnz, dtype="int32")
    data = np.zeros(nnz, dtype="float32")

    indptr_i = 1
    data_i = 0

    for cellA in range(a.shape[0]):
        for cellB in range(a.shape[0]):
            indptrA = a.indptr[cellA]
            indptrB = a.indptr[cellB]
            s = 0
            while indptrA < a.indptr[cellA+1] and indptrB < a.indptr[cellB+1]:
                if a.indices[indptrA] < a.indices[indptrB]:
                    s += square( a.data[indptrA] )
                    indptrA += 1
                elif a.indices[indptrA] > a.indices[indptrB]:
                    s += square( a.data[indptrB] )
                    indptrB += 1
                else:  # a.indices[indptrA] == a.indices[indptrB]
                    s += square( a.data[indptrA] - a.data[indptrB] )
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
        indptr[indptr_i] = data_i
        indptr_i += 1
    return(data, indices, indptr)
