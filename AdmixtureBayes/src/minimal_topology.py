from Rtree_to_covariance_matrix import get_populations
from tree_plotting import plot_node_structure_as_directed_graph
from Treemix_to_AdmixtureBayes import Node

def tree_to_node_combinations(tree):
    return get_populations(tree)

def combine_nodes(node_structure, new_node, seen_sets):
    candidate=new_node.name
    seen=[]
    for lists_of_fixed_size in seen_sets[::-1]:
        for attached_branch in lists_of_fixed_size:
            if( attached_branch.issubset(candidate) and 
               ((not attached_branch.issubset(seen)) or (not node_structure[attached_branch].has_parent()))):
                seen.extend(list(attached_branch))
                new_node.add_child(node_structure[attached_branch])
                node_structure[attached_branch].add_parent(new_node)
    return node_structure
    
def node_combination_to_node_structure(node_combinations):
    length_sorted={}
    for node_combination in node_combinations:
        leaves=frozenset(node_combination.split('.'))
        k=len(leaves)
        if k in length_sorted:
            length_sorted[k].append(leaves)
        else:
            length_sorted[k]=[leaves]
    length_sorted_list=[length_sorted.get(k,[]) for k in range(1,max(length_sorted.keys())+1)]
    #length_sorted_list is of the form [[[A],[B],[C]],[[A,B],[B,C]],...,[[A,B,C]]]
    node_structure={}
    for leaf_node in length_sorted_list[0]:
        node_structure[leaf_node]=Node(leaf_node)
    added_sets=[length_sorted_list[0]]
    for lists_of_fixed_size in length_sorted_list[1:]:
        for branch_set in lists_of_fixed_size:
            new_node=Node(branch_set)
            combine_nodes(node_structure, new_node, added_sets)
            node_structure[branch_set]=new_node
        added_sets.append(lists_of_fixed_size)
    return node_structure

def plot_node_structure(node_structure, prefix='tmp'):
    plot_node_structure_as_directed_graph(node_structure, drawing_name=prefix+'.png')
    
    