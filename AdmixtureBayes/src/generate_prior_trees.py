from numpy.random import choice
from Rtree_operations import _rename_root
from tree_plotting import pretty_print
from copy import deepcopy
from numpy import argsort

def set_outgoing_branch(node, parent_name, branch, length):
    node[branch]=parent_name
    node[branch+3]=length
    return node

factor_to_index_number={'parent_key':0, 
                        'parent_key2':1,
                        'admixture_proportion':2,
                        'branch_length':3,
                        'branch_length2':4,
                        'child_key':5,
                        'child_key2':6}

def create_node(**kwargs):
    return update_node([None]*7, **kwargs)

def update_node(node, **kwargs):
    for factor, value in kwargs.items():
        node[factor_to_index_number[factor]]=value
    return node

def generate_admix_topology(size, admixes, leaf_nodes=None):
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
    #print 'deciding allowance of:', chosen_indexes, no_totally_free
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
    
def _pair_everyhting_up_nicely(indexes, no_totally_free, halfly_frees, no_admixes, node_name):
    parents=[]
    types=[]
    seen_parents={}
    for n,index in enumerate(indexes):
        type_of_index= _classify_type(index, no_totally_free*2, len(halfly_frees), no_admixes)
        if type_of_index=='free':
            if (index-1 in seen_parents and index%2==1 and index<no_totally_free*2): 
                parents.append(seen_parents[index-1])
                types.append('second_coalescence')
            elif (index+1 in seen_parents and index%2==0 and index+1<no_totally_free*2):
                parents.append(seen_parents[index+1])
                types.append('second_coalescence')
            else:
                parent_key=node_name()
                parents.append(parent_key)
                types.append('first_coalescence')
                seen_parents[index]=parent_key
        elif type_of_index=='admix':
            parents.append(node_name(admixture=True))
            types.append('admixture')
        else:
            which_suitor=index-no_totally_free*2
            parents.append(halfly_frees[which_suitor])
            types.append('second_coalescence')
    return parents, types
            

def simulate_generation(no_totally_free, halfly_frees, no_admixes, lineages, tree, node_name):
    no_halfly_frees=len(halfly_frees)
    new_lineages=[]
    #print 'lineages', lineages
    #pretty_print(tree)
    #print 'no_totally_free', no_totally_free
    #print 'no_halfly_frees',no_halfly_frees
    #print 'no_admixes', no_admixes
    indexes=choice(no_totally_free*2+no_halfly_frees+no_admixes, size = len(lineages), replace=False )
    while not _allowed_generation(sorted(indexes), no_totally_free*2):
        indexes=choice(no_totally_free*2+no_halfly_frees+no_admixes, size = len(lineages), replace=False)
    #indexes=[1, 4, 6, 8, 9, 11]
    parent_keys, types= _pair_everyhting_up_nicely(indexes, no_totally_free, halfly_frees, no_admixes, node_name)
    #print zip(parent_keys, types)
    new_lineages=[]
    for (key,branch), parent_key, typ in zip(lineages, parent_keys, types):
        if typ=='first_coalescence':
            tree[key]=set_outgoing_branch(tree[key], parent_key, branch, 0.14)
            no_totally_free-=1
            halfly_frees.append(parent_key)
            tree[parent_key]=create_node(child_key=key)
        if typ=='second_coalescence':
            tree[key]=set_outgoing_branch(tree[key], parent_key, branch, 0.13)
            halfly_frees.remove(parent_key)
            tree[parent_key]=update_node(tree[parent_key], child_key2=key)
            new_lineages.append((parent_key,0))
        if typ=='admixture':
            no_admixes-=1
            tree[key]=set_outgoing_branch(tree[key], parent_key, branch, 0.12)
            tree[parent_key]=create_node(admixture_proportion=0.6, child_key=key)
            new_lineages.append((parent_key,0))
            new_lineages.append((parent_key,1))
    return new_lineages, tree, no_totally_free, halfly_frees, no_admixes
                    
def _classify_type(index, n_frees, n_halfs, n_admixs):
    #print index,n_frees,n_halfs,n_admixs
    if index<n_frees:
        return 'free'
    elif index<n_halfs+n_frees:
        return 'half'
    return 'admix'

def _has_partner(index, indexes):
    #print index,indexes
    if len(indexes)>0:
        print indexes[-1]
    if len(indexes)>0 and indexes[-1]==index-1 and indexes[-1]%2==0:
        return True
    return False
    
    
if __name__=='__main__':
    print _classify_type(12, 12, 0, 1)
    print _allowed_generation([1, 4, 6, 8, 9, 11],12)
    print _allowed_generation([2,3], 5)
    print _allowed_generation([1,2], 5)
    print _allowed_generation([2,3,5], 5)
    
    from tree_plotting import plot_graph, pretty_print
    ak=generate(6, 1)
    pretty_print(ak)
    plot_graph(ak)
    
    
    
    
    
    