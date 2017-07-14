from Rtree_operations import pretty_string, create_trivial_tree, get_specific_branch_lengths, update_specific_branch_lengths
from copy import deepcopy
from numpy.random import normal
from Rtree_to_coefficient_matrix import get_orthogonal_branch_space
from numpy.linalg import norm

    
def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in dic.items():
        l_determined[i]=e
    return l_determined

def get_added_branch_pieces(org, param):
    return org.dot(normal(scale=param, size=org.shape[1]))

def rescale(tree, sigma=0.01, pks={}):
    pks['rescale_constrained_adap_param']=sigma
    new_tree=deepcopy(tree)
    orgs, bi,_= get_orthogonal_branch_space(new_tree)
    #print norm(orgs, axis=0)
    branches=reverse_dic_to_list(bi)
    branch_pieces= get_added_branch_pieces(orgs, sigma)
    new_tree= update_specific_branch_lengths(new_tree, branches, branch_pieces, add=True)
    if new_tree is None:
        return tree,1,0 #rejecting by setting backward jump probability to 0.
    return new_tree ,1,1

class rescale_constrained_class(object):
    new_nodes=0
    proposal_name='rescale_constrained'
    adaption=True
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
    admixture_change=0
    reverse='rescale_constrained'
    
    def __call__(self,*args, **kwargs):
        return rescale(*args, **kwargs)


if __name__=='__main__':
    from tree_plotting import plot_graph
    from Rcatalogue_of_trees import tree_on_the_border2_with_children
    from Rtree_operations import pretty_string
    new_tree,_,_=rescale(tree_on_the_border2_with_children)
    
    print 'old_tree', pretty_string(tree_on_the_border2_with_children)
    print 'new_tree', pretty_string(new_tree)
    