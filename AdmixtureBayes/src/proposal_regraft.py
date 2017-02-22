from numpy.random import choice
from tree_operations import get_index_of_nonroot_branches, tree_prune
from copy import deepcopy

def regraft_non_root_branch(tree):
    
    index_of_nonroot_branches=get_index_of_nonroot_branches(tree)
    index_of_removed_branch, index_of_source_branch=choice(index_of_nonroot_branches,2)#choosing the non root_tree
    
    print index_of_nonroot_branches, index_of_removed_branch, index_of_source_branch
    #get_number_of_destination_branches(tree, index_of_removed_branch)
    
    remove_code=tree[index_of_removed_branch][1]
    add_code=tree[index_of_source_branch][1]
    
    if add_code in remove_code:
        #this is not possible
        return
    
    new_tree=deepcopy(tree)
    new_tree=tree_prune(tree, remove_code)
    
    
    
if __name__=="__main__":
    from tree_operations import get_tree_flatter_list
    
    tree=get_tree_flatter_list()
    
    regraft_non_root_branch(tree)
    
    