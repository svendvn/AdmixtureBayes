from posterior import initialize_posterior
from numpy import array, set_printoptions
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from Rtree_operations import get_pruned_tree_and_add, scale_tree_copy, pretty_string
from copy import deepcopy
from Rtree_to_covariance_matrix import make_covariance
from reduce_covariance import reduce_covariance
from sympy.stats.rv_interface import covariance

def get_true_posterior(wishart_file='tmp_covariance_and_multiplier.txt', 
                       true_tree_file='tmp_scaled_true_tree.txt', 
                       p=0.5, 
                       use_skewed_distr=False, 
                       wishart_df_file='tmp_wishart_DF.txt',
                       outgroup='out'):
    true_tree, nodes=read_tree_file(true_tree_file)
    set_printoptions(precision=2)
    multiplier=8.57292960745
    print make_covariance(scale_tree_copy(true_tree,multiplier), node_keys=nodes)
    print reduce_covariance(make_covariance(scale_tree_copy(true_tree,multiplier), node_keys=nodes),8)
    print true_tree
    print nodes
    print outgroup
    print pretty_string(true_tree)
    x= get_pruned_tree_and_add(true_tree, outgroup)
    print pretty_string(x[0])
    print x[1]
    
    print nodes
    nodes.remove(outgroup)
    covariance, multiplier= read_wishart_file(wishart_file, nodes)
    print 'multiplier', multiplier
    #multiplier=1.0
    print covariance
    #print make_covariance(scale_tree_copy(x[0], multiplier))
    print (make_covariance(scale_tree_copy(x[0],1.0))+x[1])*multiplier
    print (make_covariance(scale_tree_copy(x[0],1.0))+x[1])*multiplier-covariance
    print (make_covariance(scale_tree_copy(x[0],1.0))+x[1])*multiplier/covariance
    wishart_df=read_wishart_df_file(wishart_df_file)
    posterior, multiplier2=initialize_posterior(covariance, M=wishart_df, p=p, use_skewed_distr=use_skewed_distr, multiplier=multiplier, nodes=nodes)
    return posterior((scale_tree_copy(x[0], multiplier2),x[1]))
    
def read_wishart_df_file(filename):
    with open(filename, 'r') as f:
        res=float(f.readline())
    return res    

def read_tree_file(filename):
    with open(filename, 'r') as f:
        lines=f.readlines()
        print lines
        nodes=lines[0].rstrip().split()
        print nodes
        tree=identifier_to_tree_clean(lines[1].rstrip(), leaves=generate_predefined_list_string(deepcopy(nodes)))
    return tree, nodes
    
def read_wishart_file(filename, nodes):
    covariance=[]
    with open(filename, 'r') as f:
        lines=f.readlines()
        nodes2=lines[0].split()
        if len(lines)==len(nodes2)+2:#in this case, there is a multiplier
            multiplier=float(lines[-1].split('=')[-1])
        else:
            multiplier=None
        for i in range(1,len(nodes2)+1):
            covariance.append(map(float,lines[i].split()[1:]))
    covariance=array(covariance)
    if nodes!=nodes2:
        mapping={val:key for key, val in enumerate(nodes2)}
        new_order=[mapping[node] for node in nodes]
        covariance=covariance[:,new_order][new_order]
    return covariance, multiplier
            
if __name__=='__main__':
    print get_true_posterior()