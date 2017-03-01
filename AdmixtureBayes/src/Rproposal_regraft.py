from copy import deepcopy
from numpy.random import choice, random, exponential
from Rtree_operations import (get_parents, is_root, get_descendants_and_rest, 
node_is_non_admixture, node_is_non_admixture, has_child_admixture, insert_children_in_tree,
remove_parent_attachment, graft, node_is_admixture)
from random import getrandbits

def _get_possible_regrafters(tree):
    res=[]
    for key in tree:
        parents=get_parents(tree[key])
        print parents
        if len([p for p in parents if p is not None])==1:
            parent=parents[0]
            #print key,(not is_root(parent)),node_is_non_admixture(tree[parent]),(not has_child_admixture(tree, parent))
            #if (not is_root(parent)) and node_is_non_admixture(tree[parent]) and (not has_child_admixture(tree, parent)):
            if node_is_non_admixture(tree[parent]) and (not has_child_admixture(tree, parent)):
                res.append(key)
    return res

def regraft_non_root(tree):
    
    possible_nodes=_get_possible_regrafters(tree)
    
    if len(possible_nodes)==0:
        return tree, 1,1,1
    
    new_tree= deepcopy(tree)
    regrafter=choice(possible_nodes, 1)[0]
    children, other= get_descendants_and_rest(tree, regrafter)
    recipient_branch= choice(other+['r'],1)[0]
    new_tree, forward_backward= regraft(new_tree, regrafter, recipient_branch)
    _, new_other =  get_descendants_and_rest(new_tree, regrafter)
    print len(other), len(new_other)
    
    return new_tree, 1.0/len(other), 1.0/len(new_other),1

def regraft(tree, remove_key, add_to_branch):
    tree, remove_distrub, remove_val=remove_parent_attachment(tree, remove_key)
    if add_to_branch=='r':
        insertion_spot=exponential()
    else:
        insertion_spot=random()
        if node_is_admixture(tree[add_to_branch]):
            which_branch=choice(2,1)[0]
        which_branch=0
    tree=graft(tree, remove_key, add_to_branch, insertion_spot, getrandbits(128), which_branch)
    return tree
    

if __name__=='__main__':
    
    import Rtree_operations
    
    tr=insert_children_in_tree(Rtree_operations.tree_on_the_border2)
    print tr
    print _get_possible_regrafters(tr)
    print regraft_non_root(tr)