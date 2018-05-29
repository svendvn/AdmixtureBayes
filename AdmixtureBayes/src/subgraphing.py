from Rtree_to_covariance_matrix import leave_node, _thin_out_dic, Population, _merge_pops, _full_node, _add_to_waiting, follow_branch, leave_node

def leave_node(key, node, population, sub_tree):
    if node_is_non_admixture(node): 
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, target_nodes=target_nodes, child_key=key)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch_length=node[3], admixture_proportion=node[2], target_nodes=target_nodes, child_key=key, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch_length=node[4], admixture_proportion=node[0], target_nodes=target_nodes, child_key=key, dependent='none')]
        
def follow_branch(parent_key, branch_length, admixture_proportion):
    
    
def follow_branch(parent_key, branch_length, population, target_nodes, child_key, dependent="none"):
    if population
    return parent_key, population, dependent

def get_subgraph(tree, subgraph_keys):
    node_keys=tree.keys()
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=zip(node_keys,pops)
    waiting_nodes={}
    taken_nodes=[]
    sub_tree={}
    pop_strings=[]
    while True:
        for key,pop in ready_nodes:
            #pop_strings.append(pop.get_population_string(min_w))
            upds=leave_node(key, tree[key], pop, covmat)
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
            big_pop=ready_nodes[0][1]
            pop_strings.append(big_pop.get_population_string(min_w))
            break

    return sorted(list(set(pop_strings)))