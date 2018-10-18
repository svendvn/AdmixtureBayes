import numpy as np
from scipy.stats import wishart, multivariate_normal, uniform, binom
from generate_prior_covariance import generate_covariance
from reduce_covariance import reduce_covariance, Areduce
from df_estimators import variance_mean_based, likelihood_mean_based
from numpy.random import choice
from itertools import product
from covariance_estimator import initor
from Rtree_operations import scale_tree, get_leaf_keys
from construct_covariance_choices import ms_simulate_wrapper
from tree_to_data import get_xs_and_ns_from_treemix_file
from construct_filter_choices import make_filter 

def simulate_ps(N, Sigma, rho):
    res=[]
    n=Sigma.shape[0]
    Sigmai=np.linalg.inv(Sigma)
    Sigma12=np.ones((n,n))*rho
    res.append(multivariate_normal.rvs(mean=[0]*n, cov=Sigma))
    conditional_variance= Sigma-Sigma12.dot(Sigmai).dot(Sigma12)
    for i in range(1, N):
        conditional_mean=Sigma12.dot(Sigmai).dot(res[i-1])
        res.append(multivariate_normal.rvs(mean=conditional_mean, cov=conditional_variance))
    return np.array(res)

def decentralize(Pdiff, no_individuals):
    no_snps,no_individuals =Pdiff.shape
    p0s=uniform.rvs(size=no_snps)
    P=(Pdiff.T+p0s).T
    if no_individuals>0:
        P=discretize(P, no_individuals)*1.0/no_individuals
        p0s=discretize(p0s, no_individuals)*1.0/no_individuals
    return P, p0s

def draw_bootstrap_sample(P_mat, no_samples, optional_p0=None, blocksize=100):
    N,n= P_mat.shape
    #print P_mat.shape
    samples=[]
    optional_samples=[]
    blocksize=min(N,blocksize)
    for i in range(no_samples):
        indices=[range(i*blocksize,(i+1)*blocksize) for i in choice(N/blocksize,N/blocksize)]
        indices=[a for b in indices for a in b]
        samples.append(P_mat[indices,:])
        if optional_p0 is not None:
            optional_samples.append([optional_p0[ind] for ind in indices])
    if optional_samples:
        return samples, optional_samples
    return samples,None

def calc_covs(bPs, b_p0s, collapsing):
    
    if collapsing=='nothing':
        return [np.cov(cov.T-p0s) for cov, p0s in zip(bPs,b_p0s)]
    if collapsing=='Areduce':
        return [Areduce(np.cov(cov.T-np.mean(cov, axis=1))) for cov in bPs]
    if collapsing=='Rreduce':
        return [reduce_covariance(np.cov(cov.T-np.mean(cov, axis=1)),0) for cov in bPs]
    
    

def discretize(P_mat,n):
    return binom.rvs(n, p= np.clip(P_mat,0,1))

def simulate_P_from_tree(tree, no_individuals_per_population, no_snps):
    filter=make_filter(filter_type='none')
    filename_gz=ms_simulate_wrapper(tree,
                                         sample_per_pop=no_individuals_per_population,
                                         nreps=no_snps//200,
                                         theta=0.4,
                                         sites=500000,
                                         recomb_rate=1.0,
                                         full_nodes=get_leaf_keys(tree),
                                         final_pop_size=100.0,
                                         ms_file='ms.txt',
                                         treemix_file='treemix.txt',
                                         time_adjust=False)
    x,n,_=get_xs_and_ns_from_treemix_file(filename_gz, filter)
    p=x/n
    return p.T, None
    

def get_two_estimates(no_leaves, no_snps, no_individuals_per_population,
                      rho=0.0,
                      no_admixtures=None, 
                      no_bootstrap_samples=100,
                      scale_factor=0.02, 
                      blocksize=100,
                      reps=1,
                      collapsing=['nothing', 'Areduce','Rreduce'], 
                      summaries=[np.mean, np.std],
                      sim_method=['normal','ms']):
    lik_based=[]
    var_based=[]
    sim_method=initor(sim_method)
    for _ in range(reps):
        if sim_method=='normal':
            Sigma=generate_covariance(no_leaves, scale_metod='None')*scale_factor
            Ps=simulate_ps(no_snps, Sigma, rho)
            #print Ps
            #print Ps
            Ps, p0s=decentralize(Ps, no_individuals_per_population)
        elif sim_method=='ms':
            Sigma, tree=generate_covariance(no_leaves, scale_metod='None', return_tree=True)
            tree=scale_tree(tree, scale_factor)
            Ps,p0s=simulate_P_from_tree(tree, no_individuals_per_population, no_snps)
            
            
        print Ps, p0s   
        #print Ps
        #print Ps
        bPs,b_p0s= draw_bootstrap_sample(Ps, no_bootstrap_samples, optional_p0=p0s, blocksize=blocksize)
        print 'drew bootstrap samples'
        #print bPs
        covs=np.array(calc_covs(bPs, b_p0s, collapsing))
        print 'calculated covs'
        #print covs
        #print covs
        lik_based.append(likelihood_mean_based(covs))
        print 'calculated likelihood based tickets'
        print lik_based[-1]
        var_based.append(variance_mean_based(covs))
        print 'calculated variance based'
        
        print var_based[-1]
    lik_based=[lik for lik in lik_based if lik is not None]
    lik_summa=[summary(lik_based) for summary in summaries]
    var_summa=[summary(var_based) for summary in summaries]
    return lik_summa, var_summa

def make_row(varying_arguments, fixed_arguments, kwargs, output,i ):
    res_str=str(i)+' '
    for arg in varying_arguments+fixed_arguments:
        res_str+=str(kwargs[arg])+' '
    for r in output:
        for s in r:
            res_str+=str(s)+' '
    return res_str[:-1]
    

def make_header(varying_arguments, fixed_arguments, summary_list, output_types=['lik','var']):
    res='iteration '+' '.join(varying_arguments+fixed_arguments)
    summ_names=list(zip(*summary_list)[0])
    for output_type in output_types:
        res+=' '+' '.join([output_type+'_'+summ_name for summ_name in summ_names])
    return res

def make_grid_data_set(filename, 
                       no_leavess, 
                       no_snpss, 
                       no_individuals_per_populationss, 
                       rhos,
                       collapsing,
                       blocksizes,
                       scale_factors,
                       repss,
                       summary_list,
                       sim_methods):
    fixed_arguments=[]
    list_of_fixed_arguments=[]
    varying_arguments=[]
    list_of_lists_of_arguments=[]
    keys=['no_leaves', 'no_snps', 'no_individuals_per_population','rho', 'collapsing', 'blocksize', 'scale_factor','reps','sim_method']
    params=[no_leavess, no_snpss, no_individuals_per_populationss, rhos, collapsing, blocksizes, scale_factors, repss, sim_methods]
    no_grid_points=0
    for key,param in zip(keys, params):
        if len(param)==1:
            fixed_arguments.append(key)
            list_of_fixed_arguments.append(param[0])
            no_grid_points=max(no_grid_points, 1)
        else:
            no_grid_points=max(no_grid_points, len(param))
            varying_arguments.append(key)
            list_of_lists_of_arguments.append(param)
    with open(filename, 'w') as f:
        f.write(make_header(varying_arguments, fixed_arguments, summary_list)+'\n')
        for i in range(no_grid_points):
            print i, 'of', str(no_grid_points)+'...'
            kwargs={}
            for name, vals in zip(varying_arguments, list_of_lists_of_arguments):
                kwargs[name]=vals[i]
            kwargs.update({arg:val for arg,val in zip(fixed_arguments, list_of_fixed_arguments)})
            print kwargs
            output=get_two_estimates(**kwargs)
            row=make_row(varying_arguments, fixed_arguments, kwargs, output, i)
            f.write(row+'\n')

def grid_maker(*args):
    res=[]
    for i in product(*args):
        res.append(list(i))
    return map(list, zip(*res))
    
    
if __name__=='__main__':
    no_ind, no_leaves, collapsing, no_snpss= grid_maker([2,5,10,100], [5,10], ['Areduce','Rreduce'], [2000,20000])
    
    make_grid_data_set(filename='dfres.txt', 
                       no_leavess=no_leaves, 
                       no_snpss=no_snpss, 
                       no_individuals_per_populationss=no_ind, 
                       rhos=[0.0], 
                       collapsing=collapsing, 
                       blocksizes=[100], 
                       scale_factors=[0.01], 
                       repss=[1], 
                       summary_list=[('mean',np.mean), ('std',np.std)],
                       sim_methods=['ms'])
    
    
    