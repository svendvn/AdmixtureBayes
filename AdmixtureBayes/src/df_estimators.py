from scipy.stats import wishart
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt


def I_cant_believe_I_have_to_write_this_function_myself(function, lower_limit):
    old_x=lower_limit
    new_x=lower_limit*2
    old_y=function(old_x)
    new_y=function(new_x)
    sgn=1
    step_size=lower_limit
    for i in range(20):
        while new_y<old_y:
            old_x=new_x
            new_x+=sgn*step_size
            old_y=new_y
            new_y=function(new_x)
        sgn=-sgn
        step_size*=0.5
        old_x=new_x
        new_x=old_x+sgn*step_size
        old_y=new_y
        new_y=function(new_x)
    return new_x
        

def likelihood_mean_based(sample_of_matrices):
    #print sample_of_matrices
    mean_wishart=  np.mean(sample_of_matrices, axis=0)
    #print [np.linalg.det(s) for s in sample_of_matrices]
    #print 'mean_wishart', mean_wishart
    #print np.linalg.det(mean_wishart)
    #print np.linalg.matrix_rank(mean_wishart)
    
    dets=[ np.linalg.det(mat) for mat in sample_of_matrices]
    
    if (not all(det>0 for det in dets)) or np.linalg.det(mean_wishart)==0:
        print 'RETURNING NONE'
        return None
    
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
        df=df_l
        val=np.linalg.norm(var_wishart-var_rom_mean_wishart/df)**2
        #print df, '::', np.log(val)
        return np.log(val)
    

    
    
    rval=I_cant_believe_I_have_to_write_this_function_myself(penalty_function, r)
    #rval=minimize(penalty_function, 5000, bounds=[(r,None)], method='BFGS', options={'gtol':1e-15}).x[0]
    #print rval
    #print var_rom_mean_wishart*rval
    return rval
    


if __name__=='__main__':
    Sigma=np.identity(4)
    df=1043
    xs=[wishart.rvs(scale=Sigma/df, df=df) for _ in xrange(100)]
    
    print variance_mean_based(xs,df)
    print likelihood_mean_based(xs)
    