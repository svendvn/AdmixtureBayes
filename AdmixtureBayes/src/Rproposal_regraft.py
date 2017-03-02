from copy import deepcopy
from numpy.random import choice, random, exponential
from Rtree_operations import (get_parents, is_root, get_descendants_and_rest, 
node_is_non_admixture, node_is_non_admixture, has_child_admixture, insert_children_in_tree,
remove_parent_attachment, graft, node_is_admixture, get_real_parents)
from random import getrandbits

def _get_possible_regrafters(tree):
    res=[]
    for key in tree:
        parents=get_real_parents(tree[key])
        print parents
        if len(parents)==1:
            parent=parents[0]
            #print key,(not is_root(parent)),node_is_non_admixture(tree[parent]),(not has_child_admixture(tree, parent))
            #if (not is_root(parent)) and node_is_non_admixture(tree[parent]) and (not has_child_admixture(tree, parent)):
            if parent=='r' or node_is_non_admixture(tree[parent]):
                res.append(key)
    return res

def _get_possible_branches(tree, children, other):
    res=[]
    for oth in other:
        if oth in tree:
            res.append((oth,0))
            if tree[oth][1] is not None:
                res.append((oth,1))
    for child in children:
        if oth in tree:
            if tree[child][0] in other:
                res.append((child,0))
            elif tree[child][1] is not None and tree[child][1] in other:
                res.append((child,1))
    return res

def make_regraft(tree):
    
    possible_nodes=_get_possible_regrafters(tree)
    
    assert len(possible_nodes)>0, 'There were no regraft locations possible, which is strange because the root is regraftable and always look the same.'
    
    new_tree= deepcopy(tree)
    regrafter=choice(possible_nodes, 1)[0]
    new_tree, remove_distrub, remove_val=remove_parent_attachment(new_tree, regrafter)
    children, other= get_descendants_and_rest(tree, regrafter)
    candidates=_get_possible_branches(new_tree, children, other)+[('r',0)]
    ch= choice(len(candidates),1)[0]
    recipient_key, recipient_branch=candidates[ch]
    new_tree, forward_backward= regraft(new_tree, regrafter, recipient_key)
    _, new_other =  get_descendants_and_rest(new_tree, regrafter)
    print len(other), len(new_other)
    
    return new_tree

def regraft(tree, remove_key, add_to_branch):
    
    if add_to_branch=='r':
        insertion_spot=exponential()
    else:
        insertion_spot=random()
        if node_is_admixture(tree[add_to_branch]):
            which_branch=choice(2,1)[0]
        which_branch=0
    tree=graft(tree, remove_key, add_to_branch, insertion_spot, getrandbits(128), which_branch)
    return tree,1
    

if __name__=='__main__':
    
    import Rtree_operations
    
    tr=insert_children_in_tree(Rtree_operations.tree_on_the_border2)
    print tr
    print _get_possible_regrafters(tr)
    tree_final=make_regraft(tr)
    
    from tree_plotting import plot_graph
    plot_graph(tree_final)
    