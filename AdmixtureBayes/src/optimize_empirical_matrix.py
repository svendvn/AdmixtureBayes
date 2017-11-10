
import autograd.numpy as np
from autograd import grad

def loglik(Sigmai, ps, xs, ns):
    
    ps2=np.tanh(ps)*0.5+0.5
    p_start=ps2[0,:]
    p_rest=ps2[1:,:]
    p_diffs=p_start-p_rest
    n,N=p_diffs.shape
#     print 'p_start', p_start
#     print 'p_rest', p_rest
#     print 'p_diffs', p_diffs
#     print 'N n', N,n
    binomial_part= np.sum(np.log(ps2)*xs)+np.sum(np.log(1-ps2)*(ns-xs))
    normal_dett_part=float(N)/2.0*np.log(np.linalg.det(Sigmai))+0.5*n*np.sum(np.log(p_start*(1.0-p_start)))
    normal_exp_part=0
    for j in range(N):
        normal_exp_part+=np.dot(np.dot(p_diffs[:,j].T, Sigmai), p_diffs[:,j])*p_start[j]*(1.0-p_start[j])
    normal_exp_part*=-0.5
#     print binomial_part
#     print normal_dett_part
#     print normal_exp_part
    return binomial_part+normal_dett_part+normal_exp_part
    
    
def vector_to_matrices(y,n):
    square_matrix_part=y[:n**2]
    square_matrix=square_matrix_part.reshape(n,n)
    p_matrix_part=y[n**2:]
    N=len(p_matrix_part)/(n+1)
    p_matrix=p_matrix_part.reshape(n+1,N)
    return square_matrix, p_matrix

def vector_to_smatrix(y):
    n=int(np.sqrt(len(y)))
    return y.reshape(n,n)
    
def vector_to_matrix(y,n):
    N=len(y)/(n+1)
    return y.reshape(n+1,N)

def matrices_to_vector(matrices):
    m1,m2=matrices
    return np.concatenate((m1.flatten(),m2.flatten()))

def matrix_to_vector(matrix):
    return matrix.flatten()
    
def get_clean_log_lik(xs,ns, Sigma=None):
    n=ns.shape[0]-1
    def to_return(y):
        if Sigma is None:
            Sigmai, ps = vector_to_matrices(y,n)
        else:
            Sigmai=Sigma
            ps=vector_to_matrix(y, n)
        return loglik(Sigmai, ps, xs,ns)
    return to_return

def substitute_part_of_y(ps, sigma):
    

def sigma_opt(ps):
    ps2=np.tanh(ps)*0.5+0.5
    p_start=ps2[0,:]
    p_rest=ps2[1:,:]
    p_diffs=p_start-p_rest
    p_diffs_scaled=p_diffs/p_start
    return p_diffs_scaled.dot(p_diffs_scaled.T)/p_diffs_scaled.shape[1]

def full_maximization(xs,ns, initial_Sigma=None):
    initial_ps=xs/ns
    n=xs.shape[0]-1
    if initial_Sigma is None:
        initial_Sigma=np.identity(n)
    y=matrices_to_vector((initial_Sigma,np.arctanh(initial_ps)))
    f=get_clean_log_lik(xs,ns)
    gr=grad(f)
    for _ in xrange(10000):
        y+=gr(y)*0.05
        sigma,ps=vector_to_matrices(y, n)
        print f(y)
    return np.linalg.inv(sigma)
    
if __name__=='__main__':
    xs=np.array([[1.0,4.0,3.0,2.0],
                 [2.0,3.0,8.0,2.0],
                 [14.0,2.0,4.0,2.0]])
    ns=np.array([[10, 16, 10,10],
                 [10,10,20,11],
                 [23,43,8,10]])
    pstarts=np.array([0.5,0.1,0.5,1.0])
    ps=xs/ns
    Sigma=np.array([[200,-100],[-100,200]])
    print full_maximization(xs, ns)
    
    from sys import exit
    exit()
    print 'ps',ps
    print loglik(Sigma, ps,xs,ns)
    n=2
    y=matrices_to_vector((Sigma,np.arctanh(ps)))
    print 'matrices_to_vector((Sigma, ps),2)', y
    print vector_to_matrices(y, n)
    f=get_clean_log_lik(xs,ns)
    
    gr=grad(f)
    print f(y)
    print gr(y)
    for _ in xrange(1000):
        print f(y)
        y+=gr(y)*0.05
        sigma,ps=vector_to_matrices(y, n)
        print np.tanh(ps)*0.5+0.5
        print np.linalg.inv(sigma)
    s,ps=vector_to_matrices(y,n)
    print np.tanh(ps)*0.5+0.5
    print np.linalg.inv(sigma)
    #print np.tanh(res[1])
    