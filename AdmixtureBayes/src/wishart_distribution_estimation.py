from scipy.stats import wishart
from numpy import mean,var
from numpy.linalg import det, matrix_rank
from scipy.optimize import minimize
from load_data import read_data
import subprocess
from numpy.random import choice
from tree_to_data import treemix_to_cov



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

def estimate(sample_of_matrices):
    return var(sample_of_matrices, axis=0)

def get_partitions(lines, blocksize):
    list_of_lists=[]
    for i in range(0,len(lines)-blocksize, blocksize):
        list_of_lists.append(lines[i:(i+blocksize)])
    return list_of_lists

def bootstrap_indices(k):
    bootstrap_inds=choice(k, size=k, replace = True)
    return bootstrap_inds
    
    

def make_bootstrap_files(filename, blocksize=None, no_blocks=None, bootstrap_samples=None):
    assert (blocksize is not None) or (no_blocks is not None), 'Has to specify either block size or number of blocks'
    filenames=[]
    filename_reduced='.'.join(filename.split(".")[:-1])
    suffix=filename.split(".")[-1]
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if no_blocks is not None:
        blocksize=n/no_blocks
    line_sets=get_partitions(lines, blocksize)
    if bootstrap_samples is None:
        bootstrap_samples=len(line_sets)
    for i in range(bootstrap_samples):
        new_filename= filename_reduced+str(i)+'.'+suffix
        with open(new_filename, 'w') as g:
            g.write(first_line)
            bootstrap_inds=bootstrap_indices(len(line_sets))
            for i in bootstrap_inds:
                g.writelines(line_sets[i])
        new_filename_gz=new_filename+'.gz'
        subprocess.call(['gzip','-f', new_filename])
        filenames.append(new_filename_gz)
    return filenames, first_line.split()
                
def make_covariances(filenames, **kwargs):
    covs=[]
    for filename in filenames:
        print filename
        res=treemix_to_cov(filename,  **kwargs)
        covs.append(res*100)
    return covs

def estimate_degrees_of_freedom(filename, bootstrap_blocksize=100, no_blocks=None, no_bootstrap_samples=10, estimate_m=False, **kwargs):
    filenames, nodes=make_bootstrap_files(filename, blocksize=bootstrap_blocksize, no_blocks=no_blocks, bootstrap_samples=no_bootstrap_samples)
    print 'nodes', nodes
    print filenames
    covs=make_covariances(filenames, nodes=nodes, **kwargs)
    print covs[1]
    if estimate_m:
        return estimate(covs)
    return optimize(covs)
    



if __name__=='__main__':
    from numpy import array
    print estimate([array([[2,2],[2,2]]), array([[2,1],[0,-1]])])
    
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
