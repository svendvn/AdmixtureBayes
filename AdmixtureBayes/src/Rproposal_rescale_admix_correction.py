from Rtree_operations import update_all_admixtures, get_number_of_admixes, get_specific_branch_lengths, get_leaf_keys, update_specific_branch_lengths
from copy import deepcopy
from numpy.random import normal
from math import sqrt
from numpy.linalg.linalg import pinv
from numpy.linalg import matrix_rank
from Rtree_to_coefficient_matrix import make_coefficient_matrix
from scipy.stats import norm

from operator import mul


class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)

def rescale_admix_correction(tree, sigma=0.01, pks={}):
    k=get_number_of_admixes(tree)
    pks['rescale_admixtures_correction_adap_param']=sigma
    new_tree=deepcopy(tree)
    if k>0:
        updat=updater(sigma/sqrt(k))
        new_tree=update_all_admixtures(new_tree, updat)
        if new_tree is None:
            return tree,1,0 #rejecting by setting backward jump probability to 0.
    else:
        return new_tree,1,0.234 #not to have to deal with the admix=0 case, I return 0.234 such that the adaptive parameter is not affected by these special cases.
    new_tree, qforward, qbackward = getcorrection(tree, new_tree,sigma/2)
    if new_tree is None:
        return tree, 1,0
    return new_tree ,qforward, qbackward

def getcorrection(old_tree, new_tree,sigma):
    
    node_keys=sorted(get_leaf_keys(old_tree))
    
    A,_,bi1= make_coefficient_matrix(old_tree, node_keys=node_keys)
    B,_,bi2= make_coefficient_matrix(new_tree, node_keys=node_keys)
    
    branches=reverse_dic_to_list(bi1)
    
    x_A=get_specific_branch_lengths(old_tree, branches)
    x_B=get_specific_branch_lengths(new_tree, branches)
    
    lambd=pinv(B.dot(B.T)).dot(A-B).dot(x_A)
    
    mu_new=(B.T).dot(lambd)+x_A
    
    x_new=mu_new+norm.rvs(scale=sigma, size= len(mu_new))
    
    q_forward=reduce(mul, norm.pdf(mu_new-x_new, scale=sigma))
    #print 'x_A', x_A
    #print 'x_B', x_B
    #print 'x_new', x_new
    
    reverse_lambd=pinv(A.dot(A.T)).dot(B-A).dot(x_new)
    reverse_mu_new=(A.T).dot(reverse_lambd)+x_new
    
    #print 'matrix_rank , dimension (A)', matrix_rank(A), A.shape
    #print 'matrix_rank , dimension (B)', matrix_rank(B), B.shape
    #print 'x_reverse', reverse_mu_new
    
    q_backward=reduce(mul, norm.pdf(reverse_mu_new-x_A, scale=sigma))

    #wear the new values
    #print branches
    
    new_tree=update_specific_branch_lengths(new_tree,branches, x_new)
    
    

    return new_tree, q_forward, q_backward
    
    
def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in dic.items():
        l_determined[i]=e
    return l_determined

class rescale_admix_correction_class(object):
    new_nodes=0
    proposal_name='rescale_admix_correction'
    adaption=True
    input='tree'
    require_admixture=1
    reverse_require_admixture=1
    admixture_change=0
    reverse='rescale_admix_correction'
    
    def __call__(self, *args, **kwargs):
        return rescale_admix_correction(*args, **kwargs)


if __name__=='__main__':
    from tree_plotting import plot_graph
    from Rcatalogue_of_trees import tree_on_the_border2_with_children
    plot_graph(tree_on_the_border2_with_children)
    new_tree,f,b=rescale_admix_correction(tree_on_the_border2_with_children)
    print b/f
    plot_graph(new_tree)
    
    print 'tree', tree_on_the_border2_with_children
    print 'new_tree', new_tree
    