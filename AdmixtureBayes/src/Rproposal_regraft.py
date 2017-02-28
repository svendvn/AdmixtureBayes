from copy import deepcopy
from numpy.random import choice
from Rtree_operations import get_parents, is_root, get_descendants_and_rest, node_is_non_admixture, node_is_non_admixture, has_child_admixture, insert_children_in_tree

def _get_possible_regrafters(tree):
    res=[]
    for key in tree:
        parents=get_parents(tree[key])
        print parents
        if len([p for p in parents if p is not None])==1:
             
            parent=parents[0]
            #print key,(not is_root(parent)),node_is_non_admixture(tree[parent]),(not has_child_admixture(tree, parent))
            if (not is_root(parent)) and node_is_non_admixture(tree[parent]) and (not has_child_admixture(tree, parent)):
                res.append(key)
    return res

def regraft_non_root(tree):
    
    new_tree= deepcopy(tree)
    keys=tree.keys()
    
    possible_nodes=_get_possible_regrafters(tree)
        
    regrafter=choice(possible_nodes, 1)[0]
        
    children, other= get_descendants_and_rest(tree, regrafter)
    
    recipient_branch= choice(other,1)[0]
    
    new_tree= regraft(tree, regrafter, recipient_branch)
    
    _, new_other =  get_descendants_and_rest(new_tree, recipient_branch)
    print len(other), len(new_other)
    
    return new_tree, 1/len(other), 


if __name__=='__main__':
    
    import Rtree_operations
    
    tr=insert_children_in_tree(Rtree_operations.tree_on_the_border2)
    print tr
    print _get_possible_regrafters(tr)
    print regraft_non_root(tr)