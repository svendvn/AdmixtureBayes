from numpy.random import choice
from Rtree_operations import _rename_root
from tree_plotting import pretty_print
from copy import deepcopy

def generate(size, admixes, leaf_nodes=None):
    if leaf_nodes is None:
        leaf_nodes = [ 's'+str(i+1) for i in range(size)]
    free_admixes=admixes
    no_totally_free_coalescences=size-1+admixes
    no_halfly_free_coalescences=0
    
    ready_lineages=[(leaf_node,0) for leaf_node in leaf_nodes]
    tree={key:[None]*7 for key in leaf_nodes}
    admixture_node_counter=[]
    halfly_free_coalescences=[]
    
    node_name=_get_node_name()
    
    while True:
        ready_lineages, tree, no_totally_free_coalescences, halfly_free_coalescences, free_admixes = simulate_generation(no_totally_free_coalescences, 
                                                                                                                         halfly_free_coalescences, 
                                                                                                                         free_admixes, 
                                                                                                                         ready_lineages, 
                                                                                                                         tree, 
                                                                                                                         node_name)
        if no_totally_free_coalescences+len(halfly_free_coalescences)+free_admixes==0:
            break
        
    for key, node in tree.items():
        if node[0] is None and node[1] is None:
            del tree[key]
            tree=_rename_root(tree,key)
    return tree


def _allowed_generation(chosen_indexes, no_totally_free):
    print 'deciding allowance of:'
    print chosen_indexes, no_totally_free
    all_choosing_frees=all(chosen_index<no_totally_free for chosen_index in chosen_indexes)
    if not all_choosing_frees:
        return True
    if len(chosen_indexes)==1:
        return False
    no_doubles=0
    tmp=sorted(chosen_indexes)
    for c1,c2 in zip(tmp[:-1],tmp[1:]):
        if (c2==c1+1 and c1%2==0):
            no_doubles+=1
    return (no_doubles>0)

class _get_node_name(object):
    
    def __init__(self, letter='n', admixture_letter='a'):
        self.letter=letter
        self.admixture_letter=admixture_letter
        self.counter=0
        
    def __call__(self, admixture=False):
        self.counter+=1
        if admixture:
            return self.admixture_letter+str(self.counter)
        return self.letter+str(self.counter)

def simulate_generation(no_totally_free, halfly_frees, no_admixes, lineages, tree, node_name):
    no_halfly_frees=len(halfly_frees)
    new_lineages=[]
    print 'lineages', lineages
    pretty_print(tree)
    print 'no_totally_free', no_totally_free
    print 'no_halfly_frees',no_halfly_frees
    print 'no_admixes', no_admixes
    indexes=choice(no_totally_free*2+no_halfly_frees+no_admixes, size = len(lineages), replace=False )
    while not _allowed_generation(indexes, no_totally_free*2):
        indexes=choice(no_totally_free*2+no_halfly_frees+no_admixes, size = len(lineages), replace=False)
    freezed_no_frees=no_totally_free
    freezed_no_halfly=no_halfly_frees
    freezed_no_admixes=no_admixes
    freezed_halfly_frees=deepcopy(halfly_frees)
    for (key,branch), (n,index) in zip(lineages, enumerate(indexes)):
        upstream_node=_classify_type(index, freezed_no_frees*2, freezed_no_halfly, freezed_no_admixes)
        if upstream_node=='free':
            if _has_partner(index, indexes[:n]):
                brother=lineages[n-1]
                parent=tree[brother[0]][brother[1]]
                tree[key][branch]=parent
                tree[key][branch+3]=0.14
                tree[parent][6]=key
                new_lineages.append((parent,0))
                no_halfly_frees-=1
                halfly_frees.remove(parent)
            else:
                parent_key=node_name()
                tree[parent_key]=[None,None,None,None,None,key,None]
                halfly_frees.append(parent_key)
                tree[key][branch]=parent_key
                tree[key][branch+3]=0.13
                no_totally_free-=1
                no_halfly_frees+=1
        if upstream_node=='half':
            chosen_parent=freezed_halfly_frees[index-freezed_no_frees*2]
            tree[key][branch]=chosen_parent
            tree[chosen_parent][6]=key
            tree[key][branch+3]=0.12
            new_lineages.append((chosen_parent,0))
            no_halfly_frees-=1
            halfly_frees.remove(chosen_parent)
        if upstream_node=='admix':
            parent_key=node_name(admixture=True)
            tree[key][branch]=parent_key
            tree[key][branch+3]=0.11
            tree[parent_key]=[None,None,0.5,None,None,key,None]
            new_lineages.append((parent_key,0))
            new_lineages.append((parent_key,1))
            no_admixes-=1
    assert len(halfly_frees)==no_halfly_frees
    return new_lineages, tree, no_totally_free, halfly_frees, no_admixes
            
        
def _classify_type(index, n_frees, n_halfs, n_admixs):
    print index,n_frees,n_halfs,n_admixs
    if index<n_frees:
        return 'free'
    elif index<n_halfs+n_frees:
        return 'half'
    elif index<n_halfs+n_frees+n_admixs:
        return 'admix'
    else:
        assert False, 'the simulated new number was out of range'

def _has_partner(index, indexes):
    print index,indexes
    if len(indexes)>0:
        print indexes[-1]
    if len(indexes)>0 and indexes[-1]==index-1 and indexes[-1]%2==0:
        return True
    return False
    
    
if __name__=='__main__':
    print _classify_type(12, 12, 0, 1)
    print _allowed_generation([1,3,5,7,9,11],12)
    print _allowed_generation([2,3], 5)
    print _allowed_generation([1,2], 5)
    print _allowed_generation([2,3,5], 5)
    
    from tree_plotting import plot_graph, pretty_print
    ak=generate(6, 1)
    pretty_print(ak)
    plot_graph(ak)
    
    
    
    
    
    