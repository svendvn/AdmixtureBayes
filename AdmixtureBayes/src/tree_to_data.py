from Rtree_to_covariance_matrix import make_covariance
from scipy.stats import wishart
from Rtree_operations import find_rooted_nodes, get_real_parents, pretty_string, get_no_leaves, node_is_non_admixture, node_is_admixture, node_is_leaf_node, node_is_coalescence, get_real_children_root
from tree_statistics import get_timing
import subprocess
from dill.pointers import children


def tree_to_data_perfect_model(tree, df):
    m=make_covariance(tree)
    r=m.shape[0]
    m=wishart.rvs(df=r*df-1, scale=m/(r*df))
    return m

def tree_to_ms_command(tree, sample_per_pop=20, nreps=1, 
                       theta=5.0, sites=10000, recomb_rate=None):
    if recomb_rate is None:
        rec_part=' -s '+str(sites)
    else:
        rec_part=' -r '+str(recomb_rate)+ ' '+str(sites)
    n=get_no_leaves(tree)
    callstring='ms '+str(sample_per_pop*n)+' '+str(nreps)+' -t '+str(theta)+rec_part
    callstring+=' -I '+' '.join([str(sample_per_pop) for _ in xrange(n)])
    times=get_timing(tree)
    tree=extend_branch_lengths(tree,times)
    return callstring
    
def extend_branch_lengths(tree, times):
    for key, node in tree.items():
        for n,parent_key in enumerate(get_real_parents(node)):
            pseudo_time=times[parent_key]-times[key]
            tree[key][3+n]=(tree[key][3+n], pseudo_time)
    return tree

def construct_ej_en_es_string(tree, times, no_leaves):
    s_times=sorted(zip(times.values(),times.keys()))
    dic_of_lineages={(key,0):(n+1) for n,key in enumerate(s_times[:no_leaves])}
    population_count=len(dic_of_lineages)
    res_string=''
    for time,key in s_times:
        if key=='r':
            
            break
        node=tree[key]
        if node_is_leaf_node(node):
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_coalescence(node):
            i,j=get_affected_populations(tree, dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)
            dic_of_lineages[(key,0)]=j
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_admixture(node):
            
        
def get_affected_populations(tree, dic_of_lineages, children_branches):
    return dic_of_lineages[children_branches[0]], dic_of_lineages[children_branches[1]]
    
    

def trace_from_root(tree, init_freq):
    (child_key1, child_branch1, branch_length1),(child_key1, child_branch1, branch_length1)=find_rooted_nodes(tree)

if __name__=='__main__':
    from Rcatalogue_of_trees import *
    
    print tree_to_ms_command(tree_good)