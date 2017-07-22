from scipy.stats import wishart
from numpy import mean
from scipy.optimize import minimize
from load_data import read_data


def optimize(sample_of_matrices):
    
    mean_wishart= mean(sample_of_matrices, axis=0)
    r=mean_wishart.shape[0]
    
    def joint_density(df_l):
        df=df_l[0]
        val=-sum((wishart.logpdf(x, df=r*df-1, scale=mean_wishart/(r*df))  for x in sample_of_matrices))
        return val
    
    return minimize(joint_density, r, bounds=[(r,None)]).x[0]

def make_bootstrap_files(filename, blocksize=None, no_blocks=None):
    assert (blocksize is not None) or (no_blocks is not None), 'Has to specify either block size or number of blocks'
    filenames=[]
    filename_reduced=filename.split(".")[:-1]
    suffix=filename.split(".")[-1]
    with open(filename, 'r') as f:
        first_line=f.readline()
        lines=f.readlines()
    n=len(lines)
    if blocksize is not None:
        k=n/blocksize
    elif no_blocks is not None:
        blocksize=n/no_blocks
        k=no_blocks
    for i in range(k):
        new_filename= filename_reduced+str(i)+'.'+suffix
        with open(new_filename, 'w') as g:
            g.writeline(first_line)
            g.writelines(lines[:(blocksize*i)])
            g.writelines(lines[blocksize*(i+1):])
        filenames.append(new_filename)
    return filenames, first_line.split()
                
def make_covariances(filenames, **kwargs):
    covs=[]
    for filename in filenames:
        covs.append(read_data(filename, **kwargs))
    return covs

def estimate_degrees_of_freedom(filename, blocksize=1000, no_blocks=None, **kwargs):
    filenames, nodes=make_bootstrap_files(filename, blocksize=blocksize, no_blocks=no_blocks)
    covs=make_covariances(filenames, nodes=nodes, **kwargs)
    return optimize(covs)
    



if __name__=='__main__':
    from numpy import identity, ones
    matrix=identity(5)*0.9+ones((5,5))/10
    r=matrix.shape[0]
    df=100
    xs= [wishart.rvs(df=r*df-1, scale=matrix/(r*df)) for _ in xrange(100)]
    print optimize(xs)
