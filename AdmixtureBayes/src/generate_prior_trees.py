from numpy.random import choice

def generate(size, admixes, leaf_nodes=None):
    if leaf_nodes is None:
        leaf_nodes = [ 's'+str(i+1) for i in range(size)]
    free_admixes=admixes
    no_totally_free_coalescences=size-1+admixes
    no_halfly_free_coalescences=0
    
    ready_lineages=[(leaf_node,0) for leaf_node in leaf_nodes]
    tree={}
    admixture_node_counter=[]
    coalescence_node_counter=[]
    
    node_name=_get_node_name()
    
    while True:
        ready_


def _allowed_generation(chosen_indexes, no_totally_free):
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
        self.letter=n
        self.admixture_letter=admixture_letter
        self.counter=0
        
    def __call__(self, admixture=False):
        self.counter+=1
        if admixture:
            return self.admixture_letter+str(self.counter)
        return self.letter+str(self.counter)

def simulate_generation(no_totally_free, halfly_frees, no_admixes, lineages, tree, node_name):
    no_halfly_frees=len(halfly_frees)
    indexes=choice(no_totally_free+no_halfly_free+no_admixes, size = len(lineages))
    while not _allowed_generation(indexes, no_totally_free):
        indexes=choice(no_totally_free+no_halfly_free+no_admixes, size = len(lineages))
    for (key,branch), (n,choice) in zip(lineages, enumerate(indexes)):
        upstream_node=_classify_type(choice, no_totally_free, no_halfly_free, no_admixes)
        if upstram_node=='free':
            if _has_partner(choice, indexes[:n]):
                
            else:
                
            
        elif upstream_node==''
            
        
def _classify_type(choice, n_frees, n_halfs, n_admixs):
    if choice<n_frees:
        return 'free'
    elif choice<n_halfs:
        return 'half'
    return 'admix'

def _has_partner(index, indexes):
    if indexes and indexes[-1]==index-1 and indexes[-1]%2==0:
        return True
    return False
    
    
if __name__=='__main__':
    print _allowed_generation([2,3], 5)
    print _allowed_generation([1,2], 5)
    print _allowed_generation([2,3,5], 5)
    
    
    
    
    