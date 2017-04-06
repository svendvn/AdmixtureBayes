from scipy.stats import norm
from summary import Summary

class trivial_proposal(object):
    
    def __init__(self):
        pass
    
    def __call__(self, x, pks={}):
        if isinstance(x, float):
            return x+norm.rvs(),1,1,1,1,1
        elif isinstance(x, list) and isinstance(x[0], float):
            return [y+norm.rvs() for y in x],1,1
        
    
    
class Trivial_Summary(Summary):
    
    def __init__(self):
        super(Trivial_Summary,self).__init__('trivial_summary')

    def __call__(self, **kwargs):
        x=kwargs['old_tree']
        if isinstance(x,list):
            return x[0]
        return x
    
    def summary_of_phylogeny(self, tree):
        return norm.rvs()
    
    def make_histogram(self,x,a=None):
        '''
        we expect a completely normal distribution, so we use quantiles from the normal distribution.
        '''
        
        bins=[norm.ppf(q) for q in [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99]]
        super(Trivial_Summary, self).make_histogram(x=x,a=a, bins=bins)