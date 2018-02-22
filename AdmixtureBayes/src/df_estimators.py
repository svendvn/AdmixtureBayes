from scipy.stats import wishart
import numpy as np
from scipy.optimize import minimize

def likelihood_mean_based(sample_of_matrices):
    #print sample_of_matrices
    mean_wishart=  np.mean(sample_of_matrices, axis=0)
    #print [np.linalg.det(s) for s in sample_of_matrices]
    #print 'mean_wishart', mean_wishart
    #print np.linalg.det(mean_wishart)
    #print np.linalg.matrix_rank(mean_wishart)
    r=mean_wishart.shape[0]
    
    def joint_density(df_l, verbose=False):
        df=df_l[0]
        val=-sum((wishart.logpdf(x, df=df, scale=mean_wishart/df)  for x in sample_of_matrices))
        if verbose:
            print df, ':', val
        return val
    
    return minimize(joint_density, r, bounds=[(r,None)]).x[0]


def variance_mean_based(sample_of_matrices, divisor=None):
    mean_wishart=  np.mean(sample_of_matrices, axis=0)
    var_wishart= np.var(sample_of_matrices, axis=0)
    r=mean_wishart.shape[0]
    var_rom_mean_wishart=np.square(mean_wishart)+np.outer(np.diag(mean_wishart),np.diag(mean_wishart))
    #print var_wishart
    #print var_rom_mean_wishart
    def penalty_function(df_l):
        df=df_l[0]
        val=np.linalg.norm(var_wishart-var_rom_mean_wishart/df)
        #print df, '::', np.log(val)
        return np.log(val)
    rval=minimize(penalty_function, 5000, bounds=[(r,None)], method='BFGS').x[0]
    #print rval
    #print var_rom_mean_wishart*rval
    return rval
    


if __name__=='__main__':
    Sigma=np.identity(4)
    df=1043
    xs=[wishart.rvs(scale=Sigma/df, df=df) for _ in xrange(100)]
    
    print variance_mean_based(xs,df)
    print likelihood_mean_based(xs)
    