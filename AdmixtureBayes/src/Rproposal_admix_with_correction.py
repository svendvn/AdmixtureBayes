from Rproposal_admix import addadmix, deladmix
from Rtree_to_coefficient_matrix import make_coefficient_matrix
from Rtree_operations import get_leaf_keys, get_specific_branch_lengths, update_specific_branch_lengths
from numpy.linalg import pinv
from numpy.random import normal
from numpy import array

def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in dic.items():
        l_determined[i]=e
    return l_determined

def float_equal(x,y):
    return float((x-y)**2)<1e-5

def add_random_noise(vector, param=0.01):
    return array(vector) + normal(scale=param, size=org.shape[1])

def addmix_with_correction(tree,new_node_names=None,pks={}, fixed_sink_source=None, new_branch_length=None, new_to_root_length=None):
    
    added_tree,forward, backward =addadmix(tree, 
                                           new_node_names=new_node_names,
                                           pks=pks,
                                           fixed_sink_source=fixed_sink_source,
                                           new_branch_length=new_branch_length,
                                           new_to_root_length=new_to_root_length,
                                           check_opposite=False,
                                           preserve_root_distance=False)
    
    node_keys=sorted(get_leaf_keys(tree))
    
    A,_,bi1= make_coefficient_matrix(tree, node_keys=node_keys)
    B,_,bi2= make_coefficient_matrix(added_tree, node_keys=node_keys)
    
    x_A=get_specific_branch_lengths(tree, reverse_dic_to_list(bi1))
    x_B=get_specific_branch_lengths(added_tree, reverse_dic_to_list(bi2))
    
    Binverse=pinv(B)
    Ainverse=pinv(A)
    
    tilde_x_B=Binverse.dot(A.dot(x_A))
    #random_tilde_x_B=add_random_noise(tilde_x_B)
    tilde_x_A=Ainverse.dot(B.dot(tilde_x_B))
#     
#     if all((float_equal(x,y) for x,y in zip(x_A, tilde_x_A))):
#         added_tree=update_specific_branch_lengths(added_tree, bi2, tilde_x_B)
#         if added_tree is None:
#             return tree, 1,0
#         else:
#             return tree, forward, backward
#     else:
#         new_x_B=
    
    print array(x_B)
    print array(tilde_x_B)
    print array(x_A)
    print array(tilde_x_A)
    
    print B.dot(x_B)
    print A.dot(x_A)
    print B.dot(tilde_x_B)
    print A.dot(tilde_x_A)
    
    tilde2_x_A=Ainverse.dot(B.dot(x_B))
    tilde2_x_B=Binverse.dot(A.dot(tilde2_x_A))
    
    print B.dot(x_B)
    print A.dot(x_A)
    print B.dot(tilde2_x_B)
    print A.dot(tilde2_x_A)
    
    
    t=5
    #new_tree= update_specific_branch_lengths(new_tree, branches, branch_pieces, add=True)
    
    
if __name__=='__main__':
    from generate_prior_trees import generate_phylogeny
    
    tree= generate_phylogeny(4,1)
    
    addmix_with_correction(tree)