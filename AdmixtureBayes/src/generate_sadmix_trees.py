from Rtree_to_covariance_matrix import get_populations
from generate_prior_trees import generate_phylogeny, simulate_number_of_admixture_events
from copy import deepcopy
from Rtree_operations import node_is_admixture, get_number_of_admixes, pretty_string, get_number_of_leaves
from Rproposal_admix import deladmix, addadmix
from Rtree_to_coefficient_matrix import get_rank

def generate_sadmix_tree(no_leaves, no_sadmixes=None, nodes=None, p=0.5, starting_admixes=0):
    '''
    a sadmix is a Signficant ADMIXture. 
    This function generates trees with this either by rejection-accept(if starting_admixes==no_sadmixes) or by bottom up approach (starting_admixes==0).
    starting_admixes can also be somewhere in between
    '''
    
    if no_sadmixes is None:
        no_sadmixes=simulate_number_of_admixture_events(p)
    assert starting_admixes<=no_sadmixes, 'not implemented to start with more admixes than the final sadmixes'
    tree=accept_reject_generation(no_leaves, starting_admixes, nodes)
    tree=add_sadmixes(tree, no_sadmixes)
    return tree
    
def add_sadmixes(tree, final_no_sadmixes):
    k=get_number_of_admixes(tree)
    n=get_number_of_leaves(tree)
    maxrank=n*(n+1)/2
    #print pretty_string(tree)
    for i in range(k,final_no_sadmixes):
        pops=get_rank(tree)
        assert pops<maxrank, 'Admixture event number '+str(i+1)+' cant be added because the model is already maxed out'
        names=['sad'+str(i)+'a','sad'+str(i)+'b']
        candidate_tree,_,_=addadmix(tree,new_node_names=names, preserve_root_distance=False)
        candidate_pops=get_rank(candidate_tree)
        #print 'cand_res', candidate_pops, pops
        while candidate_pops<=pops:
            #print 'rejected addition'
            candidate_tree,_,_=addadmix(tree,new_node_names=names, preserve_root_distance=False)
            candidate_pops=get_rank(candidate_tree)
            #print 'cand_res', candidate_pops, pops
        tree=candidate_tree
        #print '----------'
        print pretty_string(tree)
        
    return tree


def accept_reject_generation(no_leaves, no_sadmixes, nodes=None):
    tree=generate_phylogeny(no_leaves, no_sadmixes, leaf_nodes=nodes)
    while not admixes_are_sadmixes(tree):
        tree=generate_phylogeny(no_leaves, no_sadmixes, leaf_nodes=nodes)
    return tree


def admix_is_sadmix(tree, branch, reference_pop):
    ntree,f,b=deladmix(tree,fixed_remove=branch, preserve_root_distance=False)
    pops_n=get_populations(ntree)
    return reference_pop!=pops_n
    
def admixes_are_sadmixes(tree):
    pops=None
    for key,node in tree.items():
        if node_is_admixture(node):
            if pops is None:
                pops=get_populations(tree)
            sadmix_bool=admix_is_sadmix(tree, (key,1), pops)
            if not sadmix_bool:
                return False
    return True
            
            
if __name__=='__main__':
    from Rtree_operations import pretty_string
    t=generate_sadmix_tree(10,7)
    print pretty_string(t)