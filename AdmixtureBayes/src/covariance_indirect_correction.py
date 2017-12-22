import numpy as np
from covariance_simulation import Simulator, estimate_Sigma_wrapper
from copy import deepcopy
from scipy.stats import wishart

class Sigma_proposal(object):
    
    def __init__(self, step_size=0.5,wait=6, threshold=0.4):
        self.step_size=step_size
        self.no_smallers=0
        self.wait=wait
        self.df=step_size
        self.threshold=threshold
        
    def is_wishart(self):
        if self.no_smallers>self.wait and (self.no_smallers%2)==((self.wait+1)%2):
            return True
        return False
    
    def is_wishart_era(self):
        return self.no_smallers>self.wait
    
    def __call__(self, emp_Sigma, implied_Sigma, Sigma):
        if self.is_wishart():
            print 'Check wishart'
            df=emp_Sigma.shape[0]/self.df
            prop_Sigma= wishart.rvs(df=df, scale=Sigma/df)
            return prop_Sigma
            prop_Sigma=deepcopy(Sigma)
            n=Sigma.shape[0]
            for i in range(n):
                for j in range(i,n):
                    if np.random.random()<self.threshold:
                        prop_Sigma[i,j]+=(emp_Sigma[i,j]-implied_Sigma[i,j])*self.step_size
                        if j!=i:
                            prop_Sigma[j,i]+=(emp_Sigma[i,j]-implied_Sigma[i,j])*self.step_size
        else:
            prop_Sigma=Sigma+(emp_Sigma-implied_Sigma)*self.step_size
        return prop_Sigma
        
    def choose_next_and_adapt(self, old_distance, new_distance, old_Sigma, new_Sigma):
        if new_distance>old_distance:
            
            if self.is_wishart():
                self.df*=0.5
            else:
                self.step_size*=0.8
                self.no_smallers+=1
            return old_distance, old_Sigma
        else:
            if self.is_wishart():
                self.no_smallers+=1
                self.df*=5.0
            else:
                self.step_size*=1.1
            return new_distance, new_Sigma
        
def status_print(i, proposal, emp_Sigma, implied_Sigma, Sigma, prop_Sigma, distances, prop_dist):
    to_print='Iteration '+str(i)+'\n'
    to_print+='Step size='+str(proposal.step_size)+'('+str(proposal.no_smallers)+")"+"\n"
    to_print+='Degrees of freedom='+str(emp_Sigma.shape[0]/proposal.df)+"\n"
    to_print+='{:<22} {:>8.5f}\n'.format('Old distance=', distances[-1])
    to_print+='{:<22} {:>8.5f}\n'.format('Proposed distance=', prop_dist)
    to_print+='Diagonal elements of emp Sigma and implied Sigma:'+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(emp_Sigma)])+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(implied_Sigma)])+'\n'
    to_print+='Diagonal elements of emp Sigma and proposed Sigma:'+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(prop_Sigma)])+'\n'
    to_print+=''.join(['{:<9.4f}'.format(s) for s in np.diag(Sigma)])+'\n'
    print to_print
    


def search_sigmas(xs,ns,no_its = 100, s=1, reduce_method = 'outgroup', method_of_weighing_alleles = 'outgroup_product', init_Sigma=None, Sim=None):
    nss=np.tile(ns,s)
    if Sim is None:
        Sim=Simulator(nss, reduce_method = reduce_method, method_of_weighing_alleles = method_of_weighing_alleles)
    no,N=xs.shape
    n=no-1
    emp_pijs=xs/ns
    emp_Sigma=estimate_Sigma_wrapper(emp_pijs, reduce_method=reduce_method, method_of_weighing_alleles=method_of_weighing_alleles)
    if init_Sigma is None:
        Sigma=emp_Sigma
    else:
        Sigma=init_Sigma
    implied_Sigma=Sim.get_implied_Sigma(Sigma)
    distances=[np.linalg.norm(emp_Sigma-implied_Sigma)]
    Proposal=Sigma_proposal()
    for i in range(no_its):
        prop_Sigma=Proposal(emp_Sigma, implied_Sigma, Sigma)
        try:
            implied_Sigma=Sim.get_implied_Sigma(prop_Sigma)
            prop_dist=np.linalg.norm(emp_Sigma-implied_Sigma)
        except np.linalg.linalg.LinAlgError:
            prop_dist=float('inf')
        
        status_print(i, Proposal, emp_Sigma, implied_Sigma, Sigma, prop_Sigma, distances, prop_dist)
        new_dist, new_Sigma= Proposal.choose_next_and_adapt(distances[-1], prop_dist, Sigma,prop_Sigma)
        distances.append(new_dist)
        Sigma=new_Sigma

    return Sigma
    

    
if __name__=='__main__':
#     print vector_to_symmetric_matrix([0.3]*6)
#     
#     from sys import exit
#     exit()
#     xs=np.array([[4.0,7.0,6.0,2.0,6.0],
#                  [2.0,3.0,8.0,2.0,4.0],
#                  [14.0,2.0,8.0,2.0,3.0],
#                  [10.0,23.0,4.0,3.0,11.0]])
#     ns=np.array([[10, 16, 10,10,40],
#                  [10,10,20,11,10],
#                  [23,43,8,10,6],
#                  [23,41,8,10,14]])
    from scipy.stats import uniform, norm, multivariate_normal, binom
    from brownian_motion_generation import simulate_xs_and_ns 


    Sigma = np.identity(3)*0.03+0.02
    Sigma[2,1] = 0
    Sigma[1,2] = 0
    N = 100000
    ns = np.ones((4,N))*2
    xs, p0s_temp, true_pijs = simulate_xs_and_ns(3, N, Sigma, ns, normal_xval=False)
    np.savetxt('tmp_xs.txt', xs)
    np.savetxt('tmp_p0strue.txt', p0s_temp)
    nss = np.tile(ns,1)
    Sim = Simulator(nss, reduce_method = 'outgroup', method_of_weighing_alleles = 'average_product')
    est_Sigma=search_sigmas(xs,ns, no_its=240, s=1, Sim=Sim, method_of_weighing_alleles='average_product')
    est_Sigma2=search_sigmas(xs,ns, no_its=240, s=1, init_Sigma=Sigma, Sim=Sim, method_of_weighing_alleles='average_product')
    Sim.save_to_file('tmp_')
    print est_Sigma
    print est_Sigma2
    print Sigma
    #print simulate_Sigma(Sigma, ns, reduce_method='outgroup', method_of_weighing_alleles='outgroup_product')
#     p0s_heu=[heuristic_p0(xs, ns, alpha=float(i)/10) for i in range(11)]
#     print p0s_temp[1:10], p0s_heu[1:10]
#     print xs
#     print 'simulated data'
#     #p0s=sim_p0s(xs[0,:], ns[0,:], 100)
#     print 'simulatedNp0'
#     p0s=np.zeros((N,1))
#     p0s[:,0]=p0s_temp
# #     p0s2=np.zeros((N,1))
# #     p0s2[:,0]=p0s_heu
#     print log_lik_p0s(p0s, xs,ns,Sigma)
#     print log_lik_p0s(p0s, xs,ns,Sigma+0.001)
#     #Sigma=np.identity(3)*0.08+0.04
#     print 'true_val'
#     print Sigma
#     print update_sigma(true_pijs[1:,:], p0s)
#     print update_sigma(np.clip(true_pijs[1:,:],0,1), p0s)
#     print em(xs, ns, initial_Sigma=None, p0s=p0s, maxiter=100)
#     for p0_row in p0s_heu:
#         p0s2=np.zeros((N,1))
#         p0s2[:,0]=p0_row
#         print em(xs, ns, initial_Sigma=None, p0s=p0s2, maxiter=100)
    #print em(xs, ns, initial_Sigma=None, p0s=p0s2, maxiter=100)
    #print optimize_p0_stochastic_gradient_descent(xs,ns, K=1, init_Sigma=Sigma, batch_SNP_size=1000, batch_sim_size=1, p0s=p0s)
    #print optimize_p0_llik(xs,ns, init_Sigma=Sigma)
#     pstarts=np.array([0.5,0.1,0.5,1.0])
#     ps=xs/ns
#     Sigma=np.array([[200,-100],[-100,200]])
#     trans=transform_allele_freqs(0.01)
#     print full_maximization(xs, ns, trans=trans)
#     print ps
#     
#     from sys import exit
#     exit()
#     print 'ps',ps
#     print loglik(Sigma, ps,xs,ns)
#     n=2
#     y=matrices_to_vector((Sigma,np.arctanh(ps)))
#     print 'matrices_to_vector((Sigma, ps),2)', y
#     print vector_to_matrices(y, n)
#     f=get_clean_log_lik(xs,ns)
#     
#     gr=grad(f)
#     print f(y)
#     print gr(y)
#     for _ in xrange(1000):
#         print f(y)
#         y+=gr(y)*0.05
#         sigma,ps=vector_to_matrices(y, n)
#         print np.tanh(ps)*0.5+0.5
#         print np.linalg.inv(sigma)
#     s,ps=vector_to_matrices(y,n)
#     print np.tanh(ps)*0.5+0.5
#     print np.linalg.inv(sigma)
#     #print np.tanh(res[1])
#     