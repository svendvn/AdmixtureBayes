import numpy as np
from copy import deepcopy
from covariance_estimator import Estimator


def update_pop(Sigma, hidden_states, i, freqs, vars, p0s):
    Sigma_11=Sigma[i,i]
    Sigma_21=np.delete(Sigma[:,i],i,0)
    #print Sigma_21
    Sigma_12=Sigma_21.T
    Sigma_22i=np.linalg.inv(np.delete(np.delete(Sigma, i,0),i,1))
    hidden_states_small=np.delete(hidden_states,i,0)
    #print hidden_states_small[:,1]
    N=hidden_states.shape[1]
    Sigma_prod=np.dot(Sigma_12, Sigma_22i)
    j=1
    #print p0s[j]+np.dot(Sigma_prod , hidden_states_small[:,j] - p0s[j]  )
    means=np.array([p0s[j]+np.dot(Sigma_prod , hidden_states_small[:,j] - p0s[j]  ) for j in range(N)])
    #print means
    var_base=Sigma_11+np.dot(Sigma_prod, Sigma_21.T)
    #print var_base
    vars_normal=np.array([var_base*p0*(1-p0) for p0 in p0s])
    means2=freqs[i,:].flatten()
    vars2=vars[i,:].flatten()
    #print means.shape, vars_normal.shape, means2.shape, vars2.shape
    return (means/vars_normal+means2/vars2)/(1.0/vars_normal+1.0/vars2)

def update_sigma(hidden_states, p0s):
    scaled=(hidden_states-p0s.T)/np.sqrt(p0s.T*(1.0-p0s.T))
    return np.dot(scaled,scaled.T)/hidden_states.shape[1]

def heuristic_p0(xs,ns, alpha=1.0):
    freqs=xs/ns
    freqs=np.clip(freqs,0.01,0.99)
    return freqs[0,:]*alpha+np.mean(freqs[1:], axis=0)*(1-alpha)

def initor(a):
    if not isinstance(a, basestring):
        return a[0]
    else:
        return a

def em(xs, ns, initial_Sigma=None, p0s=None, alpha=1.0, maxiter=100):
    freqs=xs[1:,:]/ns[1:,:]
    freqs=np.clip(freqs,0.01,0.99)
    print np.min(freqs)
    hidden_states=xs[1:,:]/ns[1:,:]
    hidden_states=np.clip(hidden_states,0.01,0.99)
    if p0s is None:
        p0s=heuristic_p0(xs, ns, alpha=alpha)
    p0s=p0s.clip(0.01,0.99)
    if initial_Sigma is None:
        Sigma=update_sigma(hidden_states, p0s)
    else:
        Sigma=initial_Sigma
    vars=freqs*(1.0-freqs)
    n=xs.shape[0]-1
    count=0
    old_Sigma=Sigma
    first_dist=0
    while count<maxiter:
        for i in range(n):
            Sigma=update_sigma(hidden_states, p0s)
            hidden_states[i,:]=update_pop(Sigma, hidden_states, i, freqs, vars, p0s).T
        dist=np.linalg.norm(old_Sigma-Sigma)
        print '{}: {:.6f}'.format(count,dist)
        old_Sigma=Sigma
        if dist<1e-5:
            break
            
        
        #print 'upd'
        #print Sigma
        #print p0s[:2]
        #if true_pij is not None:
            #print true_pij[:,:2]
        #print hidden_states[:,:2]
        #print freqs[:,:2]
        count+=1
        #print np.linalg.det(Sigma)
        
    return Sigma

class EmEstimator(Estimator):
    

    def __init__(self, 
                 maxiter=100,
                 alpha=1.0,
                 initial_Sigma_generator=None,
                 ):
        super(EmEstimator, self).__init__(reduce_also=True)
        self.alpha=alpha
        self.maxiter=maxiter
        self.initial_Sigma_generator=initial_Sigma_generator
        self.initialize_Sigma()

        
    def __call__(self, xs,ns):
        return em(xs, ns, initial_Sigma=self.initial_Sigma, alpha=self.alpha, maxiter=self.maxiter)
        


if __name__=='__main__':

#     xs=np.array([[4.0,7.0,6.0,2.0,6.0],
#                  [2.0,3.0,8.0,2.0,4.0],
#                  [14.0,2.0,8.0,2.0,3.0],
#                  [10.0,23.0,4.0,3.0,11.0]])
#     ns=np.array([[10, 16, 10,10,40],
#                  [10,10,20,11,10],
#                  [23,43,8,10,6],
#                  [23,41,8,10,14]])
    est=EmEstimator()
    Sigma=np.identity(3)*0.03+0.02
    Sigma[2,1]=0
    Sigma[1,2]=0
    N=10000
    ns=np.ones((4,N))*8
    from brownian_motion_generation import simulate_xs_and_ns
    xs,p0s_temp, true_pijs=simulate_xs_and_ns(3,N, Sigma, ns, normal_xval=False)
    print est(xs, ns)
   