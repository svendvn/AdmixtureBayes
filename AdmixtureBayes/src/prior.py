from scipy.stats import expon, geom
from Rtree_operations import get_all_branch_lengths, get_all_admixture_proportions,get_number_of_admixes
from math import factorial,log

def prior(tree, p=0.5):
    return 
    admixtures=get_all_admixture_proportions(tree)
    if not all(prop>=0 and prop <=1 for prop in admixtures):
        return -float('inf')
    branches=get_all_branch_lengths(tree)
    logsum=sum(map(expon.logpdf, branches))
    logsum+= geom.logpmf(len(admixtures), 1.0-p)
    return logsum

def topological_prior(tree):
    adm=get_number_of_admixes(tree)
    no_adm_nodes=2*adm
    no_nodes=len(tree)
    no_normal_nodes=no_nodes-no_adm_nodes
    assert no_normal_nodes%2==0, 'Wrong number of nodes in tree'
    n=(no_normal_nodes+2)/2
    pure_prior=non_admixed_prior(n)
    admixed_prior=get_admix_contribs(n, no_adm_nodes)
    return pure_prior+admixed_prior
    
def get_admix_contribs(n,k):
    return -log(factorial(2*(n+k)-2))+log(factorial(2*n-4))

def non_admixed_prior(n):#n is number of leaves
    return -log(factorial(2*n-3))+log(2**(n-2))+log(factorial(n-2))
    
if __name__=='__main__':
    from Rcatalogue_of_trees import *
    print topological_prior(tree_good)