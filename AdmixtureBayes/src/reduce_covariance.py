from numpy import insert, identity,ones

def reduce_covariance(covmat, subtracted_population_index):
    reducer=insert(identity(covmat.shape[0]-1), subtracted_population_index, -1, axis=1)
    #print reducer
    return reducer.dot(covmat).dot(reducer.T)

def Areduce(mat):
    A=identity(mat.shape[0])-1.0/mat.shape[0]*ones((mat.shape))
    return A.dot(mat).dot(A)


if __name__=='__main__':
    
    
    print reduce_covariance(identity(4),0
                            )