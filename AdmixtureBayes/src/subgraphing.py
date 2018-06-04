from Rtree_to_covariance_matrix import leave_node, _thin_out_dic, Population, _merge_pops, _full_node, _add_to_waiting
from Rtree_operations import node_is_non_admixture, get_leaf_keys, get_real_parents, get_real_children, rename_root, screen_and_prune_one_in_one_out, pretty_string
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from copy import deepcopy
from find_true_trees import get_unique_plottable_tree
from collections import Counter

def leave_node(key, node, population, target_nodes, follow_branch):
    if node_is_non_admixture(node): 
        return [follow_branch(parent_key=node[0],branch=0, population=population, target_nodes=target_nodes, child_key=key)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch=0, population=population, target_nodes=target_nodes, child_key=key, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch=1, population=new_pop, target_nodes=target_nodes, child_key=key, dependent='none')]
        
class follow_branch_class(object):
    
    def __init__(self, sub_graph_nodes):
        self.sub_graph_nodes=sub_graph_nodes
        self.seen_merging=False
        
    def __call__(self, parent_key, branch, population, target_nodes, child_key, dependent='none'):
        if self.seen_merging:
            return parent_key, population, dependent
        #print population.members, population.weights
        subset=population.subset_of_the_candidates(self.sub_graph_nodes)
        if subset=='partly':
            target_nodes.append((child_key, branch))
        elif subset=='all':
            self.seen_merging=True
        return parent_key, population, dependent    

def get_branches_to_keep(tree, subgraph_keys):
    node_keys=get_leaf_keys(tree)
    pops=[Population([1.0],[node]) for node in node_keys]
    follow_branch=follow_branch_class(subgraph_keys)
    ready_nodes=zip(node_keys,pops)
    waiting_nodes={}
    taken_nodes=[]
    target_nodes=[]
    #print tree
    while True:
        #print ready_nodes
        for key,pop in ready_nodes:
        
            #pop_strings.append(pop.get_population_string(min_w))
            upds=leave_node(key, tree[key], pop, target_nodes, follow_branch)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        #print 'waiting_nodes', waiting_nodes
        #print 'ready_nodes', ready_nodes
        #print 'taken_nodes', taken_nodes
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            #big_pop=ready_nodes[0][1]
            #pop_strings.append(big_pop.get_population_string(min_w))
            break
    #print 'finished tree'

    return target_nodes

def find_root_name(tree):
    parents_seen=set()
    for k in tree:
        ps=get_real_parents(tree[k])
        for p in ps:
            parents_seen.add(p)
    rootset=parents_seen-set(tree.keys())
    assert len(rootset)==1, 'wrong size of set'+str(rootset)
    return next(iter(rootset))

def prune_to_subtree(tree,branches_to_keep):
    sub_tree={}
    for key,b in branches_to_keep:
        sub_tree[key]=tree[key]
    #print 'finding root name'
    root_name=find_root_name(sub_tree)
    #print 'renaming root'
    sub_tree=rename_root(sub_tree, root_name)
    #print 'removing children'
    sub_tree=remove_empty_children(sub_tree)
    #print 'pruning'
    sub_tree=screen_and_prune_one_in_one_out(sub_tree)
    return sub_tree

def get_subtree(tree, subgraph_keys):
    branches_to_keep=get_branches_to_keep(tree, subgraph_keys)
    return prune_to_subtree(tree, branches_to_keep)

def remove_empty_children(tree):
    for k in tree:
        child_keys=get_real_children(tree[k])
        children_to_keep=[]
        for child_key in child_keys:
            if child_key in tree:
                children_to_keep.append(child_key)
        if len(child_keys)!=len(children_to_keep):
            for n,ch in enumerate(children_to_keep):
                tree[k][5+n]=ch
            for n in range(n+1,2):
                tree[k][5+n]=None
    return tree

def get_most_likely_subgraphs_list(strees, nodes, subgraph_keys, sort_nodes=True):
    if sort_nodes:
        nodes=sorted(nodes)
    topologies={}
    n=len(strees)
    for i,stree in enumerate(strees):
        if i%(n/10)==0:
            print float(i)/n
        tree=identifier_to_tree_clean(stree, leaves=generate_predefined_list_string(deepcopy(nodes)))
        sub_tree=get_subtree(tree, subgraph_keys)
        sub_stree=get_unique_plottable_tree(sub_tree)
        sub_topology, sbranch_lengths, sadmixture_proportions = sub_stree.split(';')
        branch_lengths=map(float, sbranch_lengths.split('-'))
        if len(sadmixture_proportions)>0:
            admixture_proportions=map(float, sadmixture_proportions.split('-'))
        else:
            admixture_proportions=[]
        if sub_topology in topologies:
            topologies[sub_topology][0].append(branch_lengths)
            topologies[sub_topology][1].append(admixture_proportions)
        else:
            topologies[sub_topology]=[[branch_lengths],[admixture_proportions]]
    return topologies

def get_and_save_most_likely_substrees(sub_strees, subgraph_keys, max_num=10, total_probability=1.0, prefix='',**not_needed):
    topologies={}
    for sub_stree in sub_strees:
        sub_topology, sbranch_lengths, sadmixture_proportions = sub_stree.split(';')
        branch_lengths=map(float, sbranch_lengths.split('-'))
        if len(sadmixture_proportions)>0:
            admixture_proportions=map(float, sadmixture_proportions.split('-'))
        else:
            admixture_proportions=[]
        if sub_topology in topologies:
            topologies[sub_topology][0].append(branch_lengths)
            topologies[sub_topology][1].append(admixture_proportions)
        else:
            topologies[sub_topology]=[[branch_lengths],[admixture_proportions]]
    save_top_subgraphs(topologies, subgraph_keys,max_num=max_num, total_probability=total_probability, prefix=prefix)
    
    
    
def thin_based_on_probs(N, top_keys, total_probs):
    if total_probs>=1:
        return [tk[0] for tk in top_keys]
    
    prob=0
    return_keys=[]
    for top_key, count in top_keys:
        return_keys.append(top_keys)
        prob+=float(count)/N
        if prob>=total_probs:
            return return_keys
    print 'The subtrees doesnt sum to (requested) total posterior probability'
    return return_keys


def save_to_file(topology, nodes, proportion, filename, param_val):
    snodes=sorted(nodes)
    with open(filename, 'w') as f:
        f.write(str(proportion[0])+'/'+ str(proportion[1])+' = '+str(float(proportion[0])/proportion[1])+'\n')
        f.write('\n')
        f.write(' '.join(snodes)+'\n')
        f.write(topology+'\n')
        f.write('\n')
        for row in param_val[0]:
            f.write(' '.join(map(str, row))+'\n')
        f.write('\n')
        for row in param_val[1]:
            if len(row)>0:
                f.write(' '.join(map(str, row))+'\n')
    
    
        

def save_top_subgraphs(topologies, nodes, max_num=10, total_probability=1.0, prefix='',**not_needed):
    freq_table={topology:len(nums[0]) for topology, nums in topologies.items()}
    freq_counter=Counter(freq_table)
    top_keys=freq_counter.most_common(max_num)
    N=float(sum(freq_table.values()))
    top_keys=thin_based_on_probs(N=N, 
                                 top_keys=top_keys, 
                                 total_probs=total_probability)
    for n,topology in enumerate(top_keys):
        save_to_file(topology, 
                     nodes, 
                     proportion=(freq_table[topology],N),
                     filename=prefix+'subgraph_'+'_'.join(nodes)+'-'+str(n)+'.txt', 
                     param_val=topologies[topology])
        
def read_subgraphing_dict(filename, types=['full','topological']):
    if isinstance(types, basestring):
        types=[types]
    res=[]
    with open(filename,'r') as f:
        for lin in f.readlines():
            elements=lin.split()
            arguments={}
            if '__' in elements:
                if 'full' not in types:
                    continue
                i=elements.index('__')
                outp=elements[i+1:]
                if len(outp)>0:
                    arguments['max_num']=int(outp[0])
                if len(outp)>1:
                    arguments['total_probability']=float(outp[1])
                if len(outp)>2:
                    arguments['prefix']=outp[2]
                elements=elements[:i]
            if '++' in elements:
                if 'topological' not in types:
                    continue
                i=elements.index('++')
                outp=elements[i+1:]
                if len(outp)>0:
                    arguments['max_num']=outp[0]
                if len(outp)>1:
                    arguments['prefix']=outp[1]
                elements=elements[:i]
            arguments['subgraph_keys']=elements
            res.append(arguments)
    return res
    
        
    


def prune_subtree(tree):
    keys=tree.keys()
    for key in keys:
        if key in tree:
            remove_one_in_one_out
if __name__=='__main__':
    from generate_prior_trees import generate_phylogeny
    from tree_plotting import plot_as_directed_graph
    from Rtree_operations import pretty_string
    tree={'n12': ['n13', None, None, 0.03029546, None, 'a1', 'n9'], 'n13': ['n14', None, None, 0.005071189, None, 'n12', 's10'], 'n10': ['n14', None, None, 0.008448681, None, 'a7', 'n8'], 'n11': ['n16', None, None, 0.015004037, None, 's2', 'a7'], 'n16': ['r', None, None, 0.003602743, None, 'n11', 'n15'], 'n14': ['n15', None, None, 0.033262905, None, 'n13', 'n10'], 'n15': ['n16', None, None, 0.004648582, None, 'n14', 'n6'], 's9': ['n3', None, None, 0.004760768, None, None, None], 's8': ['n6', None, None, 0.007454988, None, None, None], 's3': ['n8', None, None, 0.007021142, None, None, None], 's2': ['n11', None, None, 0.006707348, None, None, None], 's1': ['a1', None, None, 0.002359764, None, None, None], 's7': ['n9', None, None, 0.00869728, None, None, None], 's6': ['n3', None, None, 0.016867736, None, None, None], 's5': ['a2', None, None, 0.002292605, None, None, None], 's4': ['n5', None, None, 0.000109382, None, None, None], 's10': ['n13', None, None, 0.017452689, None, None, None], 'a1': ['a4', 'n12', 0.32, 0.008438453, 0.006708239, 's1', None], 'a2': ['n5', 'n6', 0.135, 0.12066482, 0.001733599, 's5', None], 'a4': ['a7', 'n9', 0.599, 0.035053964, 0.041518206, 'a1', None], 'a7': ['n10', 'n11', 0.401, 0.029920208, 0.009904888, 'a4', None], 'n8': ['n10', None, None, 0.009578977, None, 's3', 'n5'], 'n9': ['n12', None, None, 0.033982401, None, 's7', 'a4'], 'n3': ['r', None, None, 0.069890962, None, 's9', 's6'], 'n5': ['n8', None, None, 0.008113672, None, 's4', 'a2'], 'n6': ['n15', None, None, 0.002455554, None, 's8', 'a2']}
    print plot_as_directed_graph(tree)
    sub_tree=get_subtree(tree, ['s1','s2','s3'])
    print plot_as_directed_graph(sub_tree)
    print pretty_string(sub_tree)
    #plots=get_unique_plottable_tree(sub_tree)
    #print 'gotten unique_plottable'
    #print plots