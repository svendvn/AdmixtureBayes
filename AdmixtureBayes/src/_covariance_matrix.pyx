import numpy as np    # import the python access to numpy
cimport numpy as np   # import the Cython access to numpy
cimport cython        # needed for the boundscheck, etc.

@cython.boundscheck(True)
def update_matrix(np.ndarray[np.float64_t, ndim = 2] covmat,
    double branch_length,
    list weights,
    list indices):

    cdef int numElements = len(indices)
    
    cdef double [:,:] covmat_buffer = covmat
    
    cdef int i=0
    while i <numElements:
        i_=indices[i]
        covmat_buffer[i_,i_] += weights[i]**2*branch_length
        j=i+1
        while j < numElements:
            j_=indices[j]
            covmat_buffer[i_,j_]+= weights[i]*weights[j]*branch_length
            j+=1
        i+=1
            
    return covmat


    
    