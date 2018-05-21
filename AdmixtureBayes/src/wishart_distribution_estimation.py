from scipy.stats import wishart
from numpy import mean,var, savetxt, square, diag,loadtxt
from numpy.linalg import det, matrix_rank
from scipy.optimize import minimize
from load_data import read_data
import subprocess
from numpy.random import choice
from tree_to_data import treemix_to_cov, unzip, gzip
from covariance_estimator import initor
import warnings

from construct_covariance_choices import empirical_covariance_wrapper_directly
from pathos.multiprocessing import Pool
from df_estimators import variance_mean_based, likelihood_mean_based
import os
from dill.dill import FileNotFoundError



def optimize(sample_of_matrices):
    
    mean_wishart=  mean(sample_of_matrices, axis=0)
    print [det(s) for s in sample_of_matrices]
    print 'mean_wishart', mean_wishart
    print det(mean_wishart)
    print matrix_rank(mean_wishart)
    r=mean_wishart.shape[0]
    
    def joint_density(df_l, verbose=True):
        df=df_l[0]
        val=-sum((wishart.logpdf(x, df=df, scale=mean_wishart/df)  for x in sample_of_matrices))
        if verbose:
            print df, ':', val
        return val
    
    return minimize(joint_density, r, bounds=[(r,None)]).x[0]


def optimize2(sample_of_matrices):
    mean_wishart=  mean(sample_of_matrices, axis=0)
    var_wishart= var(sample_of_matrices, axis=0)
    
    var_rom_mean_wishart=(square(mean_wishart)+diag(mean_wishart)).T+diag(mean_wishart)
    
    return max(mean(var_wishart)/mean(var_rom_mean_wishart), mean_wishart.shape[0])
    
def estimate(sample_of_matrices):
    return var(sample_of_matrices, axis=0)

def get_partitions(lines, blocksize):
    list_of_lists=[]
    for i in range(0,len(lines)-blocksize, blocksize):
        list_of_lists.append(lines[i:(i+blocksize)])
    return list_of_lists

def bootstrap_indices(k):
    #print 'Bootstrapping is done '
    bootstrap_inds=choice(k, size=k, replace = True)
    return bootstrap_inds
    
    

def make_bootstrap_files(filename, blocksize=None, no_blocks=None, bootstrap_samples=None):
    assert (blocksize is not None) or (no_blocks is not None), 'Has to specify either block size or number of blocks'
    filenames=[]
    if filename.endswith('.gz'):
        filename=unzip(filename)
    filename_reduced=filename+'boot.'
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if no_blocks is not None:
        blocksize=n/no_blocks
    line_sets=get_partitions(lines, blocksize)
    print 'total SNPs=', n
    print 'total blocksize', blocksize
    print 'no_blocks', no_blocks
    print 'bootstrap_samples', bootstrap_samples
    print 'len(line_sets)', len(line_sets)
    if bootstrap_samples is None:
        bootstrap_samples=len(line_sets)
    for i in range(bootstrap_samples):
        new_filename= filename_reduced+str(i)
        with open(new_filename, 'w') as g:
            g.write(first_line)
            bootstrap_inds=bootstrap_indices(len(line_sets))
            for i in bootstrap_inds:
                g.writelines(line_sets[i])
        gzipped_filename=gzip(new_filename, overwrite=True)
        filenames.append(gzipped_filename)
    return filenames, first_line.split()

def combine_covs(tuple_covs, indices):
    cov_sum, scale_sum=0.0,0.0
    for i in indices:
        cov_sum+=tuple_covs[i][0]
        scale_sum+=tuple_covs[i][1]
    return cov_sum/scale_sum

def bootsrap_combine_covs(covs, cores, bootstrap_samples):
    p=Pool(cores)
    def t(empty):
        indices=bootstrap_indices(len(covs))
        return combine_covs(covs,indices)
    result_covs=map(t, range(bootstrap_samples))
    return result_covs                
                
def remove_files(filenames):
    for fil in filenames:
        os.remove(fil)
        os.remove(fil[:-3])                
                
def make_covariances(filenames, cores, **kwargs):
    covs=[]
    p=Pool(cores)
    def t(filename):
        return empirical_covariance_wrapper_directly(filename, **kwargs)
    covs=map(t, filenames)
    try:
        remove_files(filenames)
    except OSError as e:
        warnings.warn('Erasing the files did not succeed',UserWarning)
    return covs

def make_single_files(filename,blocksize, no_blocks):
    assert (blocksize is not None) or (no_blocks is not None), 'Has to specify either block size or number of blocks'
    filenames=[]
    if filename.endswith('.gz'):
        filename=unzip(filename)
    filename_reduced=filename+'boot.'
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if no_blocks is not None:
        blocksize=n/no_blocks
    line_sets=get_partitions(lines, blocksize)
    print 'total SNPs=', n
    print 'total blocksize', blocksize
    print 'no_blocks', no_blocks
    print 'len(line_sets)', len(line_sets)
    for i,lins in enumerate(line_sets):
        new_filename= filename_reduced+str(i)
        with open(new_filename, 'w') as g:
            g.write(first_line)
            g.writelines(lins)
        gzipped_filename=gzip(new_filename, overwrite=True)
        filenames.append(gzipped_filename)
    return filenames, first_line.split()

def estimate_degrees_of_freedom_scaled_fast(filename, 
                                            bootstrap_blocksize=1000,
                                            no_blocks=None,
                                            no_bootstrap_samples=10,
                                            summarization=['mle_opt','var_opt','opt'],
                                            cores=1,
                                            save_covs='',
                                            prefix='',
                                            load_bootstrapped_covariances=[],
                                            **kwargs):
    if not load_bootstrapped_covariances:
        single_files, nodes=make_single_files(filename, blocksize=bootstrap_blocksize, no_blocks=no_blocks)
        print single_files
        print nodes
        single_covs=make_covariances(single_files, cores=cores, return_also_mscale=True, **kwargs)
        covs=bootsrap_combine_covs(single_covs, cores=cores, bootstrap_samples=no_bootstrap_samples)
        if save_covs:
            for i,cov in enumerate(covs):
                filn=prefix+save_covs+str(i)+'.txt'
                savetxt(filn, cov)
    else:
        covs=[loadtxt(bcov) for bcov in load_bootstrapped_covariances]
    summarization=initor(summarization)
    if summarization=='var':
        res=estimate(covs)
    elif summarization=='mle_opt':
        res=likelihood_mean_based(covs)
    elif summarization=='var_opt':
        res=variance_mean_based(covs)
    return res, covs


def estimate_degrees_of_freedom(filename, 
                                bootstrap_blocksize=100, 
                                no_blocks=None, no_bootstrap_samples=10, 
                                summarization=['mle_opt','var_opt','var'],
                                cores=1, 
                                save_covs='',
                                prefix='',
                                load_bootstrapped_covariances=[],
                                **kwargs):
    if not load_bootstrapped_covariances:
        filenames, nodes=make_bootstrap_files(filename, blocksize=bootstrap_blocksize, no_blocks=no_blocks, bootstrap_samples=no_bootstrap_samples)
        print 'nodes', nodes
        print filenames
        covs=make_covariances(filenames, cores=cores, **kwargs)
        if save_covs:
            for i,cov in enumerate(covs):
                filn=prefix+save_covs+str(i)+'.txt'
                savetxt(filn, cov)
    else:
        covs=[loadtxt(bcov) for bcov in load_bootstrapped_covariances]
    #print covs[1]
    summarization=initor(summarization)
    if summarization=='var':
        res=estimate(covs)
    elif summarization=='mle_opt':
        res=likelihood_mean_based(covs)
    elif summarization=='var_opt':
        res=variance_mean_based(covs)
    return res, covs
    



if __name__=='__main__':
    
    nodes='Anzick Athabascan Greenlander Han Ket Koryak Malta  Saqqaq TA6 USR1 UstIshim'.split()
    
    from construct_filter_choices import make_filter
    locus_filter=make_filter()
    
    estimator_arguments=dict(reducer='Yoruba',
                         variance_correction='None',
                         method_of_weighing_alleles='average_sum',
                         arcsin_transform=False,
                         jade_cutoff=1e-5,
                         bias_c_weight='default',
                         nodes=nodes+['Yoruba'],
                         indirect_correction=False,
                         Indirect_its=1,
                         Indirect_multiplier_s=1,
                         EM_maxits=1,
                         EM_alpha=0.5,
                         no_repeats_of_cov_est=1,
                         Simulator_fixed_seed=True,
                         initial_Sigma_generator={'default':None},
                         locus_filter_on_simulated=locus_filter,
                         add_variance_correction_to_graph=False,
                         save_variance_correction=False,
                         prefix='sletmig/')
    
    
    estimate_degrees_of_freedom_scaled_fast(filename='sletmig/treemix_reduced2.gz', 
                                            bootstrap_blocksize=10000,
                                            no_bootstrap_samples=10,
                                            summarization='var_opt',
                                            cores=4,
                                            save_covs='bcov',
                                            prefix='sletmig/_',
                                            load_bootstrapped_covariances=[],
                                            est=estimator_arguments, 
                                            locus_filter=locus_filter)

    
    #f=estimate_degrees_of_freedom('sletmig/_treemix_in.txt', reduce_also=True, reducer='out', blocksize=100)
    #print 'f',f
    #from numpy import identity, ones
    #matrix=identity(5)*0.9 + ones((5,5))/10
    #r=matrix.shape[0]
    #df=100
    #xs= [wishart.rvs(df=r*df-1, scale=matrix/(r*df)) for _ in xrange(100)]
    #print optimize(xs)
    #blocksizes=[250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000]
    #blocksizes=[500]
    #y=[estimate_degrees_of_freedom('sletmig/_treemix_in.txt', bootstrap_blocksize=blocksize, reduce_also=True, reducer='out',  blocksize=1000) for blocksize in blocksizes]
    #print y
