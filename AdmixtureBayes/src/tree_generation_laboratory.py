from tree_plotting import plot_as_directed_graph, plot_graph, pretty_print
from tree_statistics import unique_identifier_and_branch_lengths, identifier_to_tree_clean
from generate_prior_trees import generate_phylogeny
from scipy.stats import geom
from Rproposal_admix import addadmix
from Rtree_to_covariance_matrix import make_covariance


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

def see_covariance_matrix(stree):
    print make_covariance(identifier_to_tree_clean(stree))
    
def print_tree(stree):
    pretty_print(identifier_to_tree_clean(stree))
        
    
    
    
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
    
    plot_big_tree('w.w.w.w.w.w.a.a.w-c.w.c.c.w.c.5.0.w.3.2-c.w.w.0.c.4.w-c.w.0.c.4-w.c.1-c.0;0.07-0.974-1.016-0.089-0.81-0.086-1.499-0.052-1.199-2.86-0.403-0.468-0.469-1.348-1.302-1.832-0.288-0.18-0.45-0.922-2.925-3.403;0.388-0.485')
    