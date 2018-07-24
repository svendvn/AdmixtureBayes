from tree_plotting import plot_as_directed_graph, plot_graph, pretty_print
from tree_statistics import unique_identifier_and_branch_lengths, identifier_to_tree_clean
from generate_prior_trees import generate_phylogeny
from scipy.stats import geom
from Rproposal_admix import addadmix
from Rtree_to_covariance_matrix import make_covariance
from Rtree_operations import add_outgroup, get_number_of_leaves, scale_tree_copy
from tree_to_data import reduce_covariance, tree_to_ms_command, call_ms_string, emp_cov_to_file
from minimal_topology import *

def plot_string_tree(stree):
    plot_graph(identifier_to_tree_clean(stree))
    
def plot_big_tree(stree):
    plot_as_directed_graph(identifier_to_tree_clean(stree))
    
def plot_minimal_topology(stree):
    tree=identifier_to_tree_clean(stree)
    node_combination=tree_to_node_combinations(tree)
    node_structure=node_combination_to_node_structure(node_combination)
    plot_node_structure(node_structure,'minimal')
    
def simulate_tree(no_leaves, no_admixes=None):
    if no_admixes is None:
        no_admixes=geom.rvs(p=0.5)-1
    tree=generate_phylogeny(no_leaves, no_admixes)
    return unique_identifier_and_branch_lengths(tree)

def save_tree(tree, filename):
    with open(filename, 'w') as f:
        f.write(tree)

def save_tree_and_nodes(tree, filename, nodes):
    with open(filename, 'w') as f:
        f.write(' '.join(nodes)+'\n')
        f.write(tree)

def load_tree(filename):
    with open(filename, 'r') as f:
        ls=f.readlines()
        stree=ls[-1].rstrip()
    return stree

def identifier_to_tree_clean_wrapper(stree):
    return identifier_to_tree_clean(stree)

#addadmix(tree,new_node_names=None,pks={}, fixed_sink_source=None, new_branch_length=None, new_to_root_length=None, check_opposite=False):
def add_random_admix(stree, *kwargs):
    tree=identifier_to_tree_clean(stree)
    ad = addadmix(tree, new_node_names=['x1','x2'], *kwargs)
    return unique_identifier_and_branch_lengths(ad[0])

def see_covariance_matrix(stree, reduce=None, factor=1.0):
    if reduce is None:
        return make_covariance(identifier_to_tree_clean(stree))*factor
    else:
        return reduce_covariance(make_covariance(identifier_to_tree_clean(stree)),0)*factor
    
def print_tree(stree):
    pretty_print(identifier_to_tree_clean(stree))

def get_number_of_zeros(row):
    return len([x for x in row if x>-1e-6 and x<1e-6])
    
def autogenerate_tree(no_leaves, no_admixtures, minimum_number_of_nonzeros=1, minimum_number_of_zeros=1):
    while True:
        tree= generate_phylogeny(no_leaves, no_admixtures)
        cov=make_covariance(tree)
        zeros=[get_number_of_zeros(row) for row in cov]
        no_non_zeros=cov.shape[0]-max(zeros)
        if no_non_zeros>=minimum_number_of_nonzeros and max(zeros)>=minimum_number_of_zeros:
            break
    tree=add_outgroup(tree, 'z', 0.234,1.96, 'Aa')
    cov=make_covariance(tree)
    print cov
    print reduce_covariance(cov, 0)
    plot_as_directed_graph(tree)
    suffix= str(no_leaves)+'_'+str(no_admixtures)+'_'+str(minimum_number_of_nonzeros)+'_'+str(minimum_number_of_zeros)
    return unique_identifier_and_branch_lengths(tree), suffix
    

    
    
if __name__=='__main__':
    #s=simulate_tree(5, 0)
    #save_tree(s, 'tree.txt')
    #print_tree(s)
    
    #designed_tree={
#                    's1':['n1',None, None, 0.1,None,None, None],
#                    's2':['n2',None,None,1,None,None,None],
#                    's3':['n1',None, None, 0.05, None, None, None],
#                    'out':['r',None, None, 1, None, None, None],
#                    'n1':['n2',None, None, 4, None, 's1', 's3'],
#                    'n2':['r',None, None, 1,None,'s2','n1' ]}
#     s=unique_identifier_and_branch_lengths(designed_tree)
#     save_tree(s, 'tree_d.txt')
    #print(s)
    
    #s=load_tree('tree.txt')
    #s=add_random_admix(s)
    #print_tree(s)
    #plot_big_tree(s)
    #save_tree(s, 'tree2.txt')
    
    s=load_tree('tree2.txt')
    s=add_random_admix(s)
    plot_big_tree(s)
    plot_minimal_topology(s)
    save_tree(s, 'tree3.txt')
    
    #s=load_tree('tree3.txt')
    #see_covariance_matrix(s)
    #plot_string_tree(s)
    #s, suffix=autogenerate_tree(5,2,3,1)
    #save_tree(s, 'tree'+'_'+suffix+'.txt')
    
    #s=load_tree('tree_10_2_4_2.txt')
    #print_tree(s)
    #plot_big_tree(s)
    #cov= get_empirical_matrix(s, 0.01, 20, 400)
    #print cov
    #see_covariance_matrix(s,0, 0.01)
    #emp_cov_to_file(cov, 'emp_cov_10_0_3_3.txt')
    #s_tree='a.a.a.w.w.w.a.w.w-a.a.w.w.w.w.w.w.c.a.8.w.a-a.w.w.w.w.c.c.c.w.6.w.5.7.w.w.w-w.w.w.w.c.4.c.c.w.7.w.6.w.w-w.w.w.c.c.3.w.4.w.w.w-c.w.w.0.c.w.4.w.w-c.c.w.1.w.0.w-c.c.w.0.1-c.w.0-c.0;1.003-0.423-0.833-0.055-0.898-0.624-0.171-1.274-1.223-0.724-0.686-0.105-0.202-0.05-0.757-0.109-1.003-0.467-1.328-1.234-0.541-1.566-0.11-0.549-0.287-0.099-0.015-1.245-0.068-0.826-0.335-0.423-0.357-2.784-2.184-0.919-0.656-1.384-0.187-2.054-1.742-1.522-0.035;0.393-0.67-0.322-0.292-0.385-0.195-0.393-0.549-0.467'
    #print_tree(s_tree)
    #plot_big_tree(s_tree)
    