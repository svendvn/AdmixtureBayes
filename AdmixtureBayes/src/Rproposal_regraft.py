from copy import deepcopy
from numpy.random import choice, random, exponential
from Rtree_operations import (get_parents, is_root, get_descendants_and_rest, 
node_is_non_admixture, has_child_admixture, insert_children_in_tree,
remove_parent_attachment, graft, node_is_admixture, get_real_parents, halfbrother_is_uncle)
from random import getrandbits
#from os import urandom

def _get_possible_regrafters(tree):
    res=[]
    for key in tree:
        parents=get_real_parents(tree[key])
        #print parents
        if len(parents)==1:
            parent=parents[0]
            #print key,(not is_root(parent)),node_is_non_admixture(tree[parent]),(not has_child_admixture(tree, parent))
            #if (not is_root(parent)) and node_is_non_admixture(tree[parent]) and (not has_child_admixture(tree, parent)):
            if parent=='r' or (node_is_non_admixture(tree[parent]) and not halfbrother_is_uncle(tree, key, parent)):
                res.append(key)
    return res

def _get_possible_branches(tree, children, other):
    '''
    returns the keys of all non-rooted nodes, that can be grafted into.
    '''
    res=[]
    for oth in other:
        if oth in tree:
            res.append((oth,0))
            if tree[oth][1] is not None:
                res.append((oth,1))
    for child in children:
        if child in tree:
            if tree[child][0] in other:
                res.append((child,0))
            elif tree[child][1] is not None and tree[child][1] in other :
                res.append((child,1))
    
    #removing closed branches, that is branches from admixture to root or from leaves to root. Normally, they would get the position of root, but because they are indespensible, they get a closed category.
    #this clumsy exception suggests, that the root should have its own key in the tree.
    res=[(r,w) for r,w in res if tree[r][3+w]!='closed_branch']
    
    return res

def make_regraft(tree, new_node=None):
    
    possible_nodes=_get_possible_regrafters(tree)
    
    assert len(possible_nodes)>0, 'There were no regraft locations possible, which is strange because the root is regraftable and always look the same.'
    
    new_tree= deepcopy(tree)
    regrafter=choice(possible_nodes, 1)[0]
    print 'regrafter', regrafter
    new_tree, remove_distrub, remove_val=remove_parent_attachment(new_tree, regrafter)
    children, other= get_descendants_and_rest(tree, regrafter)
    candidates=_get_possible_branches(new_tree, children, other)+[('r',0)]
    ch= choice(len(candidates),1)[0]
    recipient_key, recipient_branch=candidates[ch]
    print 'regrafter', regrafter
    print 'into_tree', candidates[ch]
    print 'new_tree',new_tree
    new_tree, forward_backward= regraft(new_tree, regrafter, recipient_key, new_node=new_node, which_branch=recipient_branch)
    _, new_other =  get_descendants_and_rest(new_tree, regrafter)
    #print len(other), len(new_other)
    
    
    
    return new_tree

def regraft(tree, remove_key, add_to_branch, new_node=None,which_branch=0):
    
    if add_to_branch=='r':
        insertion_spot=exponential()
    else:
        insertion_spot=random()
    if new_node is None:
        new_node=str(getrandbits(8)).strip()
    tree=graft(tree, remove_key, add_to_branch, insertion_spot, new_node, which_branch)
    return tree,1
    

if __name__=='__main__':
    
    import Rtree_operations
    from tree_plotting import plot_graph
    
    
#     before_illegal_tree={'a': ['n1', 'c', 0.5, 0.06604174100033824, 0.1, 's2', None], 'c': ['e', 'n2', 0.5, 0.1, 0.1, 'a', None], 'e': ['f', None, None, 0.05, None, 'c', 's3'], 'f': ['r', None, None, 0.02, None, 'n1', 'e'], 's3': ['e', None, None, 0.3, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n2', None, None, 0.1, None, None, None], 's4': ['n1', None, None, 0.3, None, None, None], 'n1': ['f', None, None, 0.18395825899966176, None, 'a', 's4'], 'n2': ['r', None, None, 1.9890637488986214, None, 'c', 's1']}
#     plot_graph(before_illegal_tree)
#     removed= remove_parent_attachment(before_illegal_tree, 'n1')[0]
#     print removed
#     nt=graft(removed, 'n1', 'a', 0.0001, 'hallo', 1)
#     plot_graph(nt)
    
    tree2={'a': ['n5', 'c', 0.5, 0.07182199688586655, 0.1, 's2', None], 'c': ['n2', 'r', 0.5, 0.15000000000000002, 0.15000000000000002, 'a', None], 'b': ['n4', None, None, 0.0032455783560232784, None, 'n5', 's4'], 's3': ['n5', None, None, 0.2930455087788524, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n4', None, None, 0.10695449122114763, None, None, None], 's4': ['b', None, None, 0.3, None, None, None], 'n2': ['r', None, None, 0.0596819405846814, None, 'n4', 'c'], 'n4': ['n2', None, None, 0.007072481059295323, None, 'b', 's1'], 'n5': ['b', None, None, 0.12817800311413347, None, 'a', 's3']}
    
    print _get_possible_regrafters(tree2)
    
    tr=insert_children_in_tree(Rtree_operations.tree_on_the_border2)
    tr1=deepcopy(tr)
    for i in range(10000):
        print 'before', tr
        #try:
        tr=make_regraft(tr, new_node='n'+str(i+1))
        #except Exception as e:
        #    print e
        #    plot_graph(tr, drawing_name=str(i)+'.png')
        #    tr=make_regraft(tr, new_node='n'+str(i+1))
        #    break
        print 'after', tr
    

    
    
    #plot_graph(tree_final)
    