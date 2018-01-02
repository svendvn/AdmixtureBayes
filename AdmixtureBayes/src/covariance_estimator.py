from numpy import array, mean, zeros, diag, sum, arcsin, sqrt
from reduce_covariance import reduce_covariance
#from optimize_empirical_matrix import full_maximization, transform_allele_freqs
from copy import deepcopy
import warnings

def initor(a):
    if not isinstance(a, basestring):
        return a[0]
    else:
        return a
    
class Estimator(object):
    
    def __init__(self, outgroup_name='', 
                 names=None, 
                 full_nodes=None, 
                 reduce_also=True,
                 arcsin_transform=False):
        self.arcsin_transform=arcsin_transform
        self.outgroup=''
        self.full_nodes=full_nodes
        self.names=names
        self.reduce_also=reduce_also
        #pass
    
        
    def get_reduce_index(self):
        n_outgroup=next((n for n, e in enumerate(self.full_nodes) if e==self.outgroup))
        return n_outgroup
    
    def __call__(self,xs,ns):
        pass
    

    
