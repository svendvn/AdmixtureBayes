from scipy.stats import gamma, uniform
from numpy import isfinite

def conditional_gamma_rvs(mean,shape, upper_limit):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    max_U= gamma.cdf(upper_limit, a=shape, scale=scale)
    drawn_U= uniform.rvs()*max_U
    x=gamma.ppf(drawn_U, a=shape, scale=scale)
    if False: #x<1e-8:
        print 'max_U', max_U
        print 'drawn_U', drawn_U
        print 'mu, sigma', (mean, variance)
        print 'shape,scale', (shape,scale)
        print 'upper_limit', upper_limit
        print 'x',x
        
    return gamma.ppf(drawn_U, a=shape, scale=scale)

def conditional_gamma_logpdf(value, mean, shape, upper_limit):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    logprob_conditional= gamma.logcdf(upper_limit, a=shape, scale=scale)
    if False: #not isfinite(logprob_conditional):
        print 'logprob_conditional', logprob_conditional
        print 'upperlimit', upper_limit
        print 'shape', shape
        print 'scale', scale
        print value
    return gamma.logpdf(value, a=shape, scale=scale) - logprob_conditional

def rvs(t,delta_L, shape=20):
    if delta_L<0:
        mean= -delta_L
        return t-conditional_gamma_rvs(mean, shape=shape, upper_limit=t)
    else:
        return t+clean_gamma_rvs(delta_L, shape=shape)

def logpdf(x,t,delta_L, shape=20):
    if delta_L < 0:
        mean=-delta_L
        simulated_part=t-x
        return conditional_gamma_logpdf(simulated_part, mean=mean, shape=shape, upper_limit=t)
    else:
        simulated_part=x-t
        return clean_gamma_logpdf(simulated_part, mean=delta_L, shape=shape)
        
    
    
def clean_gamma_rvs(mean, shape):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    return gamma.rvs(a=shape, scale=scale)

def clean_gamma_logpdf(value, mean, shape):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    return gamma.logpdf(value, a=shape, scale=scale)

    
def transform_to_shape_scale(mean, variance=None, shape=None):
    if shape is None:
        scale=variance/mean
        shape=mean/scale
    else:
        scale=mean/shape
    return shape,scale
    
if __name__=='__main__':
    from math import exp
    delta_L=-1.5
    t=2.05
    res=[]
    impossibles=0
    for _ in xrange(100):
        x= rvs(t,delta_L)
        lpdf= logpdf(x,t,delta_L)
        lpdf_back= logpdf(t,x,-delta_L)
        diff=-lpdf+lpdf_back
        print x, exp(lpdf), exp(lpdf_back)
        res.append(diff)
        
    from numpy import mean
    print min(res), max(res)
    print mean([min(1,exp(l)) for l in res if l<5])
    print len([l for l in res if l>5])
    