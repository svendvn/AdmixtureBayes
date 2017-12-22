import autograd.numpy as np
from autograd import grad
from copy import deepcopy
import sys
from scipy.stats import  binom, uniform, norm, beta
from scipy.stats import multivariate_normal
from pathos import multiprocessing
from math import log, isnan
from scipy.optimize import minimize


class transform_allele_freqs(object):
    
    def __init__(self, cutoff):
        self.cutoff=cutoff
        
    def __call__(self, p):
        return np.tanh(p)*(0.5-self.cutoff)+0.5
    
    def inverse(self, p):
        return np.arctanh((np.clip(p,0.4,0.6)))
    
transf=transform_allele_freqs(0.01)


def get_key(x,n, reverse=False):
    if reverse:
        return tuple( (xi,ni) for xi,ni in zip(x,n))
    else:
        return tuple( (ni-xi,ni) for xi,ni in zip(x,n))

def simmed_ps(Sigma, x,n, K):
    nm=max(n)
    p=np.clip(x/n,1/float(nm+1),float(nm-1)/float(nm+1))
    Omega=np.diag([pi*(1-pi)/ni for ni,pi in zip(n,p)])
    Omegai=np.linalg.inv(Omega)
    Omega_s=Omega[1:,1:]
    Omegai_s=Omegai[1:,1:]
    obs=np.array(p[1:])
    Sigmai=np.linalg.inv(Sigma)
    Var=np.linalg.inv(Omegai_s+Sigmai)
    print Var
    res=np.zeros((K,len(n)))
    qs=[]
    for i in range(K):
        p0=norm.rvs(p[0], scale=Omega[0,0])
        q=norm.pdf(p0, p[0], Omega[0,0])
        p0_vec=np.array([p0]*(len(x)-1))
        m=np.dot(Var, np.dot(Omegai_s, obs)+np.dot(Sigmai, p0_vec))
        p_rest=multivariate_normal.rvs(mean= m,cov=Var) 
        q*=multivariate_normal.pdf(p_rest, mean=m, cov=Var)
        res[i,:]=np.insert(p_rest,0,p0)
        qs.append(q)
    return res, qs
    

def sim_ps(xs, ns, Sigma,K):
    
    n,N=xs.shape
    dic_of_sims={}
    for col in xrange(N):
        x=xs[:,col]
        n=ns[:,col]
        print x,n
        key=get_key(x,n)
        if key in dic_of_sims:
            dic_of_sims[key]=(dic_of_sims[key][0]+1, dic_of_sims[key][1])
            continue
        rkey=get_key(x, n, reverse=True)
        if rkey in dic_of_sims:
            dic_of_sims[rkey]=(dic_of_sims[rkey][0]+1, dic_of_sims[rkey][1])
            continue
        ps, qs= simmed_ps(Sigma, x,n,K)
        dic_of_sims[key]=(1, ps,qs)
    return dic_of_sims

def sim_p0s(x0s, n0s, K):
    res=np.zeros((len(x0s),K))
    for i,(x,n) in enumerate(zip(x0s,n0s)):
        res[i,:]=beta.rvs(a=x+1, b=n-x+1, size=K)
    return res

def log_lik_p0s(p0s, xs,ns,Sigma):
    n,N = xs.shape
    res=0
    for col in range(N):
        x=xs[:,col]
        n=ns[:,col]
        p0s_row=p0s[col,:]
        res+=np.log(lik_single(p0s_row, x, n, Sigma))
    #print Sigma, res
    return res

def mvn2(cov_mat, x,mean):
    det_part=np.linalg.det(2*np.pi*cov_mat)**(-0.5)
    exp_part=np.dot(x-mean, np.dot(np.linalg.inv(cov_mat), x-mean))
    return det_part*np.exp(-0.5*exp_part)

def lik_single(p0s, x,n,Sigma):
    p=np.array(x/n)[1:]
    Omega=np.diag([pi*(1-pi)/ni for ni,pi in zip(n[1:],p)])
    res=0
    for p0 in p0s:
        p0_vec=np.array([p0]*len(p))
        res+=mvn2(Sigma*p0*(1-p0)+Omega, p, p0_vec)
    return res/len(p0s)

def get_clean_p0_llik(xs,ns,p0s, sign=1.0, SNP_from_to=[0,None],sim_from_to=[0,None]):
    if SNP_from_to[1] is None:
        SNP_from_to[1]=p0s.shape[0]
    if sim_from_to[1] is None:
        sim_from_to[1]=p0s.shape[1]
    f,t=SNP_from_to
    fs,ts=sim_from_to
    def llik(Sigma_vec):
        Sigma=vector_to_smatrix(Sigma_vec)
        return sign*log_lik_p0s(p0s[f:t,fs:ts], xs[:,f:t], ns[:,f:t], Sigma)
    return llik

def optimize_p0_stochastic_gradient_descent(xs,ns, K=100, init_Sigma=None, batch_SNP_size=1, batch_sim_size=None,p0s=None):

    if batch_sim_size is None:
        batch_sim_size=K
    x0s=xs[0,:]
    n0s=ns[0,:]
    if p0s is None:
        p0s=sim_p0s(x0s, n0s, K)
    n=xs.shape[0]-1
    no_SNPs=xs.shape[1]
    if init_Sigma is None:
        init_Sigma=(np.identity(n)-0.1)+0.1
    y=matrix_to_vector(init_Sigma)
    lik=get_clean_p0_llik(xs,ns, p0s, sign=1.0)
    glik=grad(lik)
    liks,gliks=get_partial_lik0s(no_SNPs, batch_SNP_size, K, batch_sim_size, xs, ns, p0s, sign=1.0)
    return sgd(liks, gliks, y=y, big_eval=lik, evals=10000)
    
def sgd(funs, gfuns,y, big_eval=None, evals=10000):
    step_size=0.0001
    no_funs=len(funs)
    randomized_order=np.random.permutation(no_funs)
    count=1
    while count<evals:
        i=randomized_order[count%no_funs]
        grad_val=gfuns[i](y)
        l_before=funs[i](y)
        new_y=y+grad_val*step_size
        l_after=funs[i](new_y)
        if (not np.isnan(l_after)) and l_after>l_before and check_symmetry(new_y):
            step_size*=1.0+1.0/(count+1)
            y=new_y
        else:
            step_size*=1.0-1.0/(count+1)**0.5
        print 'step_size',step_size
        if count%no_funs==0: #different cykles to prevent repitition
            randomized_order=np.random.permutation(no_funs)
            print 'finished one round'
            if big_eval is not None:
                print 'l(', str(y),')=',big_eval(y)
            else:
                print 'y=', str(y)
        count+=1
        
def check_symmetry(y):
    m=vector_to_smatrix(y)
    ans=np.array_equal(m,m.T)      
    #print ans
    return ans
    
def get_partial_lik0s(no_SNPs, batch_SNP_size, K, batch_sim_size, xs,ns,p0s, sign):
    liks=[]
    gliks=[]
    for snp in range(0, no_SNPs, batch_SNP_size):
        SNP_from_to=[snp, snp+batch_SNP_size]
        for sim in range(0, K, batch_sim_size):
            batch_from_to=[sim,sim+batch_sim_size]
            l=get_clean_p0_llik(xs, ns, p0s, sign=sign, SNP_from_to=SNP_from_to, sim_from_to=batch_from_to)
            gl=grad(l)
            liks.append(l)
            gliks.append(gl)   
    return liks,gliks
    

def optimize_p0_llik(xs,ns, K=100, init_Sigma=None):
    x0s=xs[0,:]
    n0s=ns[0,:]
    p0s=sim_p0s(x0s, n0s, K)
    n=xs.shape[0]-1
    if init_Sigma is None:
        init_Sigma=(np.identity(n)-0.1)+0.1
    y=symmetric_matrix_to_vector(init_Sigma)
    f=get_clean_p0_llik(xs,ns, p0s, sign=-1.0)
    gr=grad(f)
    return minimize(f, x0=y,  method='BFGS')
    
    step_size=0.01
    x_old=f(y)
    lik_increased=True
    last_computed_gradient=0
    index_for_last_computed_gradient=-1
    i=0
    while i<1000:
        print str(i)+'.','step_size', step_size
        z=deepcopy(y)
        #try:
        if lik_increased:
            last_computed_gradient=gr(y)
            y+=last_computed_gradient*step_size
            last_computed_gradient=i
        else:
            y+=last_computed_gradient*step_size
        #except:
        #    print 'linalg error'
        #    y=z
        #    step_size*=0.5
        #    continue
        x=f(y)
        if x_old>x or isnan(x):
            y=z
            step_size*=0.5
            lik_increased=False
        else:
            lik_increased=True
            step_size*=2
            i+=1
            x_old=x
        print x, x_old, lik_increased
        print vector_to_smatrix(y)
    
    return vector_to_smatrix(y)
    
        

    
def lik(dic_of_sims, Sigma):
    lsum=0
    for xn in dic_of_sims:
        count, ps,qs=dic_of_sims[xn]
        x,n=map(list, list(zip(*xn)))
        weights=[]
        for row in range(len(qs)):
            tsum=1
            p=ps[row,:]
            if any(pi<0 for pi in p):
                continue
            tsum*=1#np.prod(binom.pmf(x, n, p))
            tsum*=multivariate_normal.pdf(p[1:], np.array([p[0]]*(len(p)-1)), cov=Sigma)
            tsum/=qs[row]
            weights.append(tsum)
        print 'ESS=', sum(weights)**2/sum([w**2 for w in weights])
        lsum+=count*sum(weights)
    return lsum
                            
        
    
    
    


def loglik(Sigmai, ps, xs, ns, trans=transf):
    
    ps2=trans(ps)
    p_start=ps2[0,:]
    p_rest=ps2[1:,:]
    p_diffs=p_rest-p_start
    p_diffs_scaled=p_diffs/np.sqrt(p_start*(1.0-p_start))
#     np.savetxt('dump_ps.txt', ps)
#     np.savetxt('dump_ps2.txt', ps2)
#     np.savetxt('dump_xs.txt', xs)
#     np.savetxt('dump_ns.txt', ns)
#     np.savetxt('dump_pstart.txt', p_start)
#     np.savetxt('dump_pdiffs.txt', p_diffs)
#     np.savetxt('dump_Sigmai.txt', Sigmai)
    n,N=p_diffs.shape
#     print 'p_start', p_start
#     print 'p_rest', p_rest
#     print 'p_diffs', p_diffs
#     print 'N n', N,n
    binomial_part= np.sum(np.log(ps2)*xs)+np.sum(np.log(1-ps2)*(ns-xs))
    normal_dett_part=float(N)/2.0*np.log(np.linalg.det(Sigmai))-0.5*n*np.sum(np.log(p_start*(1.0-p_start)))
    normal_exp_part=0
    for j in range(N):
        normal_exp_part+=np.dot(np.dot(p_diffs_scaled[:,j].T, Sigmai), p_diffs_scaled[:,j])
    normal_exp_part*=-0.5
#     print binomial_part
#     print normal_dett_part
#     print normal_exp_part
    return binomial_part+normal_dett_part+normal_exp_part

def loglik2(Sigmai, ps, xs, ns, trans=transf):
    ps2=trans(ps)
    p_start=ps2[0,:]
    p_rest=ps2[1:,:]
    p_diffs=p_rest-p_start
    p_diffs_scaled=p_diffs/np.sqrt(p_start*(1.0-p_start))
    n,N=p_diffs.shape
    
    binomial_part=np.sum(binom.logpmf(xs, ns, ps2))
    normal_part=0
    for j in range(N):
        p0=p_start[j]
        normal_part+=multivariate_normal.logpdf(p_diffs[:,j], mean=np.array([0]*n), cov= np.linalg.inv(Sigmai)*p0*(1-p0))
    
    return binomial_part+normal_part

def loglik3(Sigma, p0s, xs,ns):
    
    ps2=trans(ps)
    p_start=ps2[0,:]
    p_rest=ps2[1:,:]
    p_diffs=p_rest-p_start
    p_diffs_scaled=p_diffs/np.sqrt(p_start*(1.0-p_start))
#     np.savetxt('dump_ps.txt', ps)
#     np.savetxt('dump_ps2.txt', ps2)
#     np.savetxt('dump_xs.txt', xs)
#     np.savetxt('dump_ns.txt', ns)
#     np.savetxt('dump_pstart.txt', p_start)
#     np.savetxt('dump_pdiffs.txt', p_diffs)
#     np.savetxt('dump_Sigmai.txt', Sigmai)
    n,N=p_diffs.shape
#     print 'p_start', p_start
#     print 'p_rest', p_rest
#     print 'p_diffs', p_diffs
#     print 'N n', N,n
    binomial_part= np.sum(np.log(ps2)*xs)+np.sum(np.log(1-ps2)*(ns-xs))
    normal_dett_part=float(N)/2.0*np.log(np.linalg.det(Sigmai))-0.5*n*np.sum(np.log(p_start*(1.0-p_start)))
    normal_exp_part=0
    for j in range(N):
        normal_exp_part+=np.dot(np.dot(p_diffs_scaled[:,j].T, Sigmai), p_diffs_scaled[:,j])
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

def vector_to_symmetric_matrix(y):
    n=(-1+np.sqrt(1+len(y)*8))/2
    print n
    n=int(n)
    res=np.zeros((n,n))
    s=0
    for i in range(n):
        for j in range(n-i):
            res[i, j]=y[s+j]
        s+=n-i
    return res+res.T-np.diag(np.diag(res))

def symmetric_matrix_to_vector(mat):
    n=mat.shape[0]
    y=[]
    for i in range(n):
        y.extend(mat[i,i:])
    return y
    
def vector_to_matrix(y,n):
    N=len(y)/(n+1)
    return y.reshape(n+1,N)



    

def update_sigma(hidden_states, p0s):
    scaled=(hidden_states-p0s.T)/np.sqrt(p0s.T*(1.0-p0s.T))
    return np.dot(scaled,scaled.T)/hidden_states.shape[1]


def matrices_to_vector(matrices):
    m1,m2=matrices
    return np.concatenate((m1.flatten(),m2.flatten()))

def matrix_to_vector(matrix):
    return matrix.flatten()
    
def get_clean_log_lik(xs,ns, Sigma=None, trans=transf, use_loglik2=False):
    n=ns.shape[0]-1
    if use_loglik2:
        logl=loglik2
    else:
        logl=loglik
    def to_return(y):
        if Sigma is None:
            Sigmai, ps = vector_to_matrices(y,n)
        else:
            Sigmai=Sigma
            ps=vector_to_matrix(y, n)
        return logl(Sigmai, ps, xs,ns, trans=trans)
    return to_return

def substitute_part_of_y(y, n, trans=transf):
    sigma,ps=vector_to_matrices(y,n)
    print 'DETERMINANT', np.linalg.det(sigma)
    print 'ps',trans(ps)
    iv=sigma_opt(ps, trans=trans)
    print 'DETERMINANT', np.linalg.det(iv)
    opt_sigma=np.linalg.inv(iv)
    return matrices_to_vector((opt_sigma, ps))



def sigma_opt(ps, trans=transf):
    ps2=trans(ps)
    p_start=ps2[0,:]
    p_rest=ps2[1:,:]
    p_diffs=p_rest-p_start
    p_diffs_scaled=p_diffs/np.sqrt(p_start*(1-p_start))
    return p_diffs_scaled.dot(p_diffs_scaled.T)/p_diffs_scaled.shape[1]

def get_one_dim_likelihood(Sigmai, xs,ns, trans=transf):
    n=len(xs)-1
    def lik(ps):
        
        ps2=trans(ps)
        p_start=ps2[0]
        p_rest=ps2[1:]
        p_diffs=p_rest-p_start
        p_diffs_scaled=p_diffs/np.sqrt(p_start*(1.0-p_start))

        binomial_part= np.sum(np.log(ps2)*xs)+np.sum(np.log(1-ps2)*(ns-xs))
        normal_dett_part=-0.5*n*np.log(p_start*(1.0-p_start))
        normal_exp_part=-0.5*np.dot(np.dot(p_diffs_scaled, Sigmai), p_diffs_scaled.T)
    #     print binomial_part
    #     print normal_dett_part
    #     print normal_exp_part
        return binomial_part+normal_dett_part+normal_exp_part
    return lik


def max_one_SNP(Sigmai, xs,ns, init_val=None,its=100, trans=transf):
    if init_val is None:
        init_val=xs/ns
    f=get_one_dim_likelihood(Sigmai, xs, ns, trans=trans)
    g=grad(f)
    step_size=0.01
    last_grad_calc=0
    last_rejected=False
    x_old=init_val
    y_old=f(x_old)
    y_first=y_old
    i=0
    tot_its=0
    while i<its and tot_its<its*2:
        if not last_rejected:
            last_grad_calc=g(x_old)
        x=x_old+last_grad_calc*step_size
        y=f(x)
        if y<=y_old:
            last_rejected=True
            step_size*=0.5
            tot_its+=1
            continue
        else:
            tot_its+=1
        last_rejected=False
        step_size*=2
        y_old=y
        x_old=x
        i+=1
        #print y_old
    #print '----------------'
        
    #print 'improvement', y_first, y_old
    return x_old

    

def full_maximization(xs,ns,initial_Sigma=None, trans=transf, it_size=20, updates=20):
    ps=trans.inverse(xs/ns)
    n=xs.shape[0]-1
    N=xs.shape[1]
    if initial_Sigma is None:
        Sigmai=np.linalg.inv((np.identity(n)-0.1)+0.1)
    else:
        Sigmai=np.linalg.inv(inital_Sigma)
    for _ in range(updates):
        dic={}
        for j in range(N):
            to_be_replaced=tuple(xs[:,j])
            if to_be_replaced in dic:
                ps[:,j]=dic[to_be_replaced]
            else:
                n=max_one_SNP(Sigmai, xs[:,j],ns[:,j], init_val=ps[:,j],its=it_size, trans=trans)
                ps[:,j]=n
                #print len(dic)
                dic[to_be_replaced]=n
        print 'DETERMINANT', np.linalg.det(Sigmai)
        iv=sigma_opt(ps, trans=trans)
        print 'ps',ps
        print 'iv',iv
        print 'DETERMINANT', np.linalg.det(iv)
        Sigmai=np.linalg.inv(iv)
    return np.linalg.inv(Sigmai)
    
    

def full_max2(xs,ns, initial_Sigma=None, trans=transf):
    initial_ps=xs/ns
    n=xs.shape[0]-1
    if initial_Sigma is None:
        initial_Sigma=(np.identity(n)-0.1)+0.1
    y=matrices_to_vector((initial_Sigma,trans.inverse(initial_ps)))
    f=get_clean_log_lik(xs,ns, trans=trans)
    f2=get_clean_log_lik(xs,ns, trans=trans, use_loglik2=True)
    gr=grad(f)
    step_size=0.01
    x_old=f(y)
    lik_not_increase=False
    last_computed_gradient=0
    index_for_last_computed_gradient=-1
    i=0
    while i<1000:
        print str(i)+'.','step_size', step_size
        z=deepcopy(y)
        #try:
        if i>last_computed_gradient:
            last_computed_gradient=gr(y)
            y+=gr(y)*step_size
            last_computed_gradient=i
        else:
            y+=last_computed_gradient*step_size
        #except:
        #    print 'linalg error'
        #    y=z
        #    step_size*=0.5
        #    continue
        x=f(y)
        if x_old>x:
            y=z
            step_size*=0.5
            continue
        else:
            step_size*=2
            i+=1
        if i%2==0:
            z=deepcopy(y)
            
            if lik_not_increase:
                pass
            
            #try:            
            y=substitute_part_of_y(y,n, trans=trans)
            #except:
            #    print 'ERROR', sys.exc_info()[0]
            #    y=z
            #    continue
            x_new=f(y)
            
            #print x,x_new
            if x_new<x:
                y=z
                lik_not_increase=True
                print 'update did not increase likelihood'
                print 'x_new',x_new
                print 'x',x
                continue
            else:
                x=x_new
            print 's',
        print x, x_old
        print x#, f2(y)
        x_old=x
        sigma,ps=vector_to_matrices(y, n)
        print np.linalg.inv(sigma)
        
        
    sigma,ps=vector_to_matrices(y, n)
    print trans(ps)
        
            
    return np.linalg.inv(sigma)

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

    Sigma=np.identity(3)*0.03+0.02
    Sigma[2,1]=0
    Sigma[1,2]=0
    N=10000
    ns=np.ones((4,N))*8
    xs,p0s_temp, true_pijs=simulate_xs_and_ns(3,N, Sigma, ns, normal_xval=False)
    est_Sigma=search_sigmas(xs,ns, no_its=100)
    print est_Sigma
    print Sigma
    print simulate_Sigma(Sigma, ns, reduce_method='outgroup', method_of_weighing_alleles='outgroup_product')
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
