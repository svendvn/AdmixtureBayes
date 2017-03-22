from scipy.stats import expon, geom
from Rtree_operations import (get_all_branch_lengths, get_all_admixture_proportions, get_number_of_admixes, get_number_of_leaves, 
get_leaf_keys,get_destination_of_lineages, get_categories)
from math import log, factorial
from scipy.special import binom 

def prior(tree, p=0.5):
    admixtures=get_all_admixture_proportions(tree)
    if not all(prop>=0 and prop <=1 for prop in admixtures):
        return -float('inf')
    branches=get_all_branch_lengths(tree)
    if not all(branch>0 for branch in branches):
        return -float('inf')
    branch_prior=sum(map(expon.logpdf, branches))
    no_admix_prior=geom.logpmf(len(admixtures)+1, 1.0-p)
    admix_prop_prior=0
    top_prior=topological_prior(tree)
    logsum=branch_prior+no_admix_prior+admix_prop_prior+top_prior
    print 'branch_prior', branch_prior
    print 'no_admix_prior',no_admix_prior
    print 'admix_prop_prior',admix_prop_prior
    print 'top_prior', top_prior
    return logsum

def topological_prior(tree):
    no_adms=get_number_of_admixes(tree)
    no_leaves=get_number_of_leaves(tree)
    return -get_admixture_extra_tops(no_leaves,no_adms)
    
def trees_without_admixture(no_leaves):
    '''
    calculating the number of rooted binary
    '''
    return -( log(factorial(2*no_leaves-3)) +
              log(2*no_leaves-3) -
              log(factorial(2**(no_leaves-2))) -
              log(factorial(no_leaves-2)))
    
def get_admixture_extra_tops(n,k):
    if k==0:
        return 0
    else:
        return sum(2*log(2*n-2+j*3)-log(2) for j in range(k))
    
def _events_selection(sampled_unoccupied, sampled_occupied, sampled_admixtures, total_unoccupied, total_occupied, total_admixtures):
    return (log(binom(total_unoccupied,sampled_unoccupied))+
            log(binom(total_occupied,sampled_occupied))+
            log(binom(total_admixtures, sampled_admixtures))-
            log(binom(total_unoccupied+total_occupied+total_admixtures,
                      sampled_unoccupied+sampled_occupied+sampled_admixtures)))
    
def _totally_free_selection():

def tree_prior(tree):
    total_prob=0
    no_admix=get_number_of_admixes(tree)
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in leaves]
    totally_free_coalescences=coalescences_nodes+['r']
    halfly_free_coalescences=[]
    free_admixtures=admixture_nodes
    no_free_coalescences=len(totally_free_coalescences)-1+no_admix
    
    coalescences_on_hold=[]
    sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
    waiting_coalescences, awaited_coalescences, still_on_hold = _matchmake(single_coalescences, coalescences_on_hold)
    p1=_events_selection(sampled_unoccupied=2*len(sames)+len(waiting_coalescences), 
                        sampled_occupied=len(awaited_coalescences), 
                        sampled_admixtures=len(admixtures),
                        total_unoccupied=2*len(totally_free_coalescences),
                        total_occupied=len(halfly_free_coalescences),
                        total_admixtures=len(free_admixtures))
    p_doubles=totally_free_selection(len(totally_free_coalescences),
    
    return sames, single_coalescences, admixtures, waiting_coalescences,awaited_coalescences,still_on_hold,p1
    #waiting,awaited=thin_out_coalescences(single_coalescences, coalescences_on_hold)
    
def _matchmake(single_coalescences, coalescences_on_hold):
    happy_couples=[]
    continuing_singles=[]
    for key,branch in coalescences_on_hold:
        if (key,branch) in single_coalescences:
            happy_couples.append(((key,branch),single_coalescences[(key,branch)]))
            del single_coalescences[(key,branch)]
        else:
            continuing_singles.append((key,branch))
    return single_coalescences, happy_couples, continuing_singles
    
    

if __name__=='__main__':
    from math import exp
    print exp(get_admixture_extra_tops(3, 2))
    from Rcatalogue_of_trees import tree_good, tree_clean, tree_one_admixture, tree_two_admixture
    from Rtree_operations import insert_children_in_tree
    print prior(tree_good)
    print prior(insert_children_in_tree(tree_clean))
    print prior(insert_children_in_tree(tree_one_admixture))
    print prior(insert_children_in_tree(tree_two_admixture))
    
    res=tree_prior(tree_good)
    for el in res:
        print el
    