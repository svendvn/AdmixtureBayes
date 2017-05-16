from Rtree_to_covariance_matrix import make_covariance
from scipy.stats import wishart
from Rtree_operations import find_rooted_nodes, get_real_parents, pretty_string, get_no_leaves, node_is_non_admixture, node_is_admixture, node_is_leaf_node, node_is_coalescence, get_real_children_root
from tree_statistics import get_timing
import subprocess
from numpy import loadtxt, cov, array, mean, vstack
from copy import deepcopy


def tree_to_data_perfect_model(tree, df):
    m=make_covariance(tree)
    r=m.shape[0]
    m=wishart.rvs(df=r*df-1, scale=m/(r*df))
    return m

def tree_to_ms_command(rtree, sample_per_pop=100, nreps=1, 
                       theta=0.1, sites=1000000, recomb_rate=None):
    tree=deepcopy(rtree)
    if recomb_rate is None:
        rec_part=' -s '+str(sites)
    else:
        rec_part=' -r '+str(recomb_rate)+ ' '+str(sites)
    n=get_no_leaves(tree)
    callstring='ms '+str(sample_per_pop*n)+' '+str(nreps)+' -t '+ str(theta)+' ' +rec_part + ' '
    callstring+=' -I '+str(n)+' '+' '.join([str(sample_per_pop) for _ in xrange(n)])+' '
    times=get_timing(tree)
    tree=extend_branch_lengths(tree,times)
    callstring+=construct_ej_en_es_string(tree, times, n)
    return callstring
    
def extend_branch_lengths(tree, times):
    for key, node in tree.items():
        for n,parent_key in enumerate(get_real_parents(node)):
            pseudo_time=times[parent_key]-times[key]
            tree[key][3+n]=(tree[key][3+n], pseudo_time)
    return tree



def construct_ej_en_es_string(tree, times, no_leaves):
    s_times=sorted(zip(times.values(),times.keys()))
    dic_of_lineages={(key,0):(n+1) for n,(time,key) in enumerate(s_times[:no_leaves])}
    population_count=len(dic_of_lineages)
    res_string=''
    for time,key in s_times:
        if key=='r':
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
            break
        node=tree[key]
        if node_is_leaf_node(node):
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_coalescence(node):
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_admixture(node):
            population_count+=1
            i=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))[0]
            res_string+='-es '+str(time)+' '+str(i)+' '+str(1.0-node[2])+' '
            dic_of_lineages[(key,0)]=i
            dic_of_lineages[(key,1)]=population_count
            pop_size1=calculate_pop_size(node[3])
            pop_size2=calculate_pop_size(node[4])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size1)+' '
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,1)])+' '+str(pop_size2)+' '
    return res_string
            
        
def get_affected_populations(dic_of_lineages, children_branches):
    return [dic_of_lineages[children_branch] for children_branch in children_branches] 
    
def calculate_pop_size(tup):
    drift, actual=tup
    return actual/drift/16

def call_ms_string(ms_string, sequence_file):
    with open(sequence_file, 'w') as f:
        p = subprocess.Popen(ms_string, stdout=subprocess.PIPE, shell=True)
        line_number = 0
        for line in p.stdout.readlines():
            line_number += 1
            if line_number >= 7 and '//' not in line:
                #print len(line)
                print line_number
                f.write(line)
            else:
                print line
    return 0

def calculate_covariance_matrix(file='tmp.txt', samples_per_pop=20, no_pops=4):
    data=[]
    with open(file, 'r') as f:
        for r in f.readlines():
            print r[:5]
            data.append(map(int,list(r.rstrip())))
    m= array(data)
    ps=tuple([mean(m[(i*samples_per_pop):((i+1)*samples_per_pop-1), : ], axis=0) for i in xrange(no_pops)])
    p=vstack(ps)
    print p
    print cov(p)

def trace_from_root(tree, init_freq):
    (child_key1, child_branch1, branch_length1),(child_key1, child_branch1, branch_length1)=find_rooted_nodes(tree)

if __name__=='__main__':
    from Rcatalogue_of_trees import *
    from Rtree_operations import create_trivial_tree
    tree2=create_trivial_tree(3)
    print call_ms_string(tree_to_ms_command(tree2, recomb_rate=500.0), 'tmp.txt')
    calculate_covariance_matrix('tmp.txt', 100, 3)
    print make_covariance(tree2)