from tree_plotting import plot_as_directed_graph, plot_graph, pretty_print
from tree_statistics import unique_identifier_and_branch_lengths, identifier_to_tree_clean
from generate_prior_trees import generate_phylogeny
from scipy.stats import geom
from Rproposal_admix import addadmix
from Rtree_to_covariance_matrix import make_covariance
from Rtree_operations import add_outgroup, get_number_of_leaves, scale_tree_copy
from tree_to_data import reduce_covariance, tree_to_ms_command, call_ms_string, ms_to_treemix2, emp_cov_to_file

def plot_string_tree(stree):
    plot_graph(identifier_to_tree_clean(stree))
    
def plot_big_tree(stree):
    plot_as_directed_graph(identifier_to_tree_clean(stree))
    
def simulate_tree(no_leaves, no_admixes=None):
    if no_admixes is None:
        no_admixes=geom.rvs(p=0.5)-1
    tree=generate_phylogeny(no_leaves, no_admixes)
    return unique_identifier_and_branch_lengths(tree)

def save_tree(tree, filename):
    with open(filename, 'w') as f:
        f.write(tree)

def load_tree(filename):
    with open(filename, 'r') as f:
        stree=f.readline().rstrip()
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
        print make_covariance(identifier_to_tree_clean(stree))*factor
    else:
        print reduce_covariance(make_covariance(identifier_to_tree_clean(stree)),0)*factor
    
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
    #s=simulate_tree(9, 0)
    #save_tree(s, 'tree.txt')
    #plot_big_tree(s)
    
    #s=load_tree('tree.txt')
    #s=add_random_admix(s)
    #print_tree(s)
    #plot_big_tree(s)
    #save_tree(s, 'tree2.txt')
    
    #s=load_tree('tree2.txt')
    #s=add_random_admix(s)
    #plot_big_tree(s)
    #save_tree(s, 'tree3.txt')
    
    #s=load_tree('tree3.txt')
    #see_covariance_matrix(s)
    #plot_string_tree(s)
    s, suffix=autogenerate_tree(5,2,3,1)
    save_tree(s, 'tree'+'_'+suffix+'.txt')
    
    #s=load_tree('tree_10_0_3_3.txt')
    #print_tree(s)
    #plot_big_tree(s)
    #cov= get_empirical_matrix(s, 0.01, 20, 400)
    #print cov
    #see_covariance_matrix(s,0, 0.01)
    #emp_cov_to_file(cov, 'emp_cov_10_0_3_3.txt')
    #s_tree='w.c.1.c.c.4.w.w.3.w-w.c.c.1.w.2.w-c.w.0.w.w-c.w.w.0-c.w.0-c.0;0.001-0.001-0.001-0.0-0.003-0.003-0.003-0.002-0.001-0.0-0.002-0.002-0.001-0.003-0.0-0.001-0.001-0.001;'
    #print_tree(s_tree)
    #plot_big_tree(s_tree)
    