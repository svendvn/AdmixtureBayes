from Rtree_to_covariance_matrix import leave_node, _thin_out_dic, Population, _merge_pops, _full_node, _add_to_waiting
from Rtree_operations import node_is_non_admixture, get_leaf_keys, get_real_parents, rename_root, screen_and_prune_one_in_one_out

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
    while True:
        print ready_nodes
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

def create_subtree(tree,branches_to_keep):
    sub_tree={}
    for key,b in branches_to_keep:
        sub_tree[key]=tree[key]
    root_name=find_root_name(sub_tree)
    sub_tree=rename_root(sub_tree, root_name)
    sub_tree=screen_and_prune_one_in_one_out(sub_tree)
    return sub_tree


def prune_subtree(tree):
    keys=tree.keys()
    for key in keys:
        if key in tree:
            remove_one_in_one_out
if __name__=='__main__':
    from generate_prior_trees import generate_phylogeny
    from tree_plotting import plot_as_directed_graph
    from Rtree_operations import pretty_string
    tree=generate_phylogeny(10,2)
    print plot_as_directed_graph(tree)
    a=get_branches_to_keep(tree, ['s1','s2','s3'])
    print a
    sub_tree=create_subtree(tree,a)
    print pretty_string(sub_tree)
    
    plot_as_directed_graph(sub_tree, drawing_name='tt.png')
    
    