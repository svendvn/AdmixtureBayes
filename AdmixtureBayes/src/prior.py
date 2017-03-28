from scipy.stats import expon, geom
from Rtree_operations import (get_all_branch_lengths, get_all_admixture_proportions, get_number_of_admixes, get_number_of_leaves, 
get_leaf_keys,get_destination_of_lineages, get_categories, get_parent_of_branch, propagate_married, propagate_admixtures)
from math import log, factorial,exp
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
    #print 'branch_prior', branch_prior
    #print 'no_admix_prior',no_admix_prior
    #print 'admix_prop_prior',admix_prop_prior
    #print 'top_prior', top_prior
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
    #print total_admixtures, total_occupied, total_unoccupied
    #print sampled_admixtures, sampled_occupied, sampled_unoccupied
    
    logmhypergeom= (log(binom(total_unoccupied,sampled_unoccupied))+
            log(binom(total_occupied,sampled_occupied))+
            log(binom(total_admixtures, sampled_admixtures))-
            log(binom(total_unoccupied+total_occupied+total_admixtures,
                      sampled_unoccupied+sampled_occupied+sampled_admixtures)))
    
    return logmhypergeom
    
def _totally_free_selection(total_pairs,sampled_pairs, sampled_ones):
    denom=log(binom(total_pairs*2,sampled_pairs*2+sampled_ones))
    num=log(2)*sampled_ones
    num+=log(factorial(total_pairs))
    num-=sum(log(factorial(x)) for x in [total_pairs-sampled_ones-sampled_pairs,
                                         sampled_pairs,
                                         sampled_ones])
    return num-denom

def _awaited_selection(positions, taken):
    return -log(binom(positions,taken))

def _assigned_selection(total, pairs, unlabeleld_singles, labelled_singles, admixtures):
    return -(log(factorial(total))-
             log(factorial(pairs))-
             log(factorial(unlabeleld_singles))-
             log(factorial(admixtures))-
             log(2)*pairs)
    
def _get_number_of_admixture_branches(lineages):
    return sum(1 for key, branch in lineages if branch==1)
    

def tree_prior(tree):
    total_prob=0
    no_admix=get_number_of_admixes(tree)
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in leaves]
    totally_free_coalescences=coalescences_nodes+['r']
    free_admixtures=admixture_nodes
    
    coalescences_on_hold=[]
    res=0
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        #print sames, waiting_coalescences, admixtures
        res+=_get_selection_probabilities(no_sames=len(sames), 
                                          no_waiting_coalescences=len(waiting_coalescences), 
                                          no_awaited_coalescences=len(awaited_coalescences), 
                                          no_admixtures=len(admixtures), 
                                          no_totally_free_coalescences=len(totally_free_coalescences), 
                                          no_free_admixtures=len(free_admixtures), 
                                          no_ready_lineages=len(ready_lineages), 
                                          no_coalescences_on_hold=len(coalescences_on_hold),
                                          no_admixture_pairs=_get_number_of_admixture_branches(ready_lineages))
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        totally_free_coalescences, free_admixtures=_thin_out_frees(tree, totally_free_coalescences,free_admixtures, ready_lineages)
        #print 'free_admixture', free_admixtures
        #print 'totally_free_coalescences', totally_free_coalescences
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))
        #print 'ready_lineages', ready_lineages

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
        
    return res
    
    
    return sames, single_coalescences, admixtures, waiting_coalescences,awaited_coalescences,still_on_hold,p1,p2,p3,p4, exp(p1+p2+p3+p4)
    #waiting,awaited=thin_out_coalescences(single_coalescences, coalescences_on_hold)

def _thin_out_frees(tree, totally_free_coalescences,free_admixtures, lineages):
    for key, branch in lineages:
        parent_key=get_parent_of_branch(tree,key,branch)
        if parent_key in totally_free_coalescences:
            totally_free_coalescences.remove(parent_key)
        if parent_key in free_admixtures:
            free_admixtures.remove(parent_key)
    return totally_free_coalescences, free_admixtures

def _get_selection_probabilities(no_sames, no_waiting_coalescences, no_awaited_coalescences,
                                 no_admixtures, no_totally_free_coalescences,
                                 no_free_admixtures, no_ready_lineages, no_coalescences_on_hold,
                                 no_admixture_pairs):
    p1=_events_selection(sampled_unoccupied=2*no_sames+no_waiting_coalescences, 
                            sampled_occupied=no_awaited_coalescences, 
                            sampled_admixtures=no_admixtures,
                            total_unoccupied=2*no_totally_free_coalescences,
                            total_occupied=no_coalescences_on_hold,
                            total_admixtures=no_free_admixtures)
    p2=_totally_free_selection(no_totally_free_coalescences,no_sames,no_waiting_coalescences)
    total_chosen=2*no_sames+no_waiting_coalescences+no_awaited_coalescences+no_admixtures
    if no_totally_free_coalescences>=total_chosen:
        p1b=_events_selection(sampled_unoccupied=total_chosen, 
                            sampled_occupied=0, 
                            sampled_admixtures=0,
                            total_unoccupied=2*no_totally_free_coalescences,
                            total_occupied=no_coalescences_on_hold,
                            total_admixtures=no_free_admixtures)
        p2b=_totally_free_selection(no_totally_free_coalescences, 0, total_chosen)
        if p1b+p2b<0:
            p1-=log(1-exp(p1b+p2b))
    p3=_awaited_selection(no_coalescences_on_hold,no_awaited_coalescences)
    p4=_assigned_selection(no_ready_lineages, no_sames, no_waiting_coalescences, no_awaited_coalescences, no_admixtures)
    return p1+p2+p3+p4

def matchmake(single_coalescences, coalescences_on_hold):
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
    #print prior(tree_good)
    
    #print prior(insert_children_in_tree(tree_one_admixture))
    #print prior(insert_children_in_tree(tree_two_admixture))
    
    #res=tree_prior(insert_children_in_tree(tree_clean))
    #print res, exp(res)

    
    from generate_prior_trees import generate
    from tree_plotting import pretty_print
    
    tree_trouble={'s3': ['n1', None, None, 0.13, None, None, None], 's2': ['n2', None, None, 0.14, None, None, None], 's1': ['n1', None, None, 0.14, None, None, None], 'a3': ['r', 'r', 0.6, 0.14, 0.13, 'n2', None], 'n1': ['n2', None, None, 0.13, None, 's1', 's3'], 'n2': ['a3', None, None, 0.12, None, 's2', 'n1']}
    print tree_prior(tree_trouble)
    
    for i in range(200):
        tree=generate(2,1)
        print tree
        pretty_print(tree)
        print exp(tree_prior(tree)),
        
    
    
    import sys
    sys.exit()
    
    sampled_ones=[]
    print sum(log(2) for i in sampled_ones)
    
    print exp(_totally_free_selection(5,0,3))
    print exp(_totally_free_selection(5, 1, 1))
    
    tree_2e2a={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['c',None,None, 0.3,None,None,None],
      'a':['b','d', 0.5,0.2,0.1,'s2',None],
      'c':['e','f',0.5,0.1,0.1,'s4',None],
      'b':['r',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['b',None,None,0.05,None,'s1','c']}
    
    tree_3e1a=tree_good
    
    tree_4e={'s1':['d',None, None, 0.1,None,None,None],
      's2':['f',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','s2'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}
    
    tree_2e2c={'s1':['d',None, None, 0.1,None,None,None],
      's2':['e',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'s2','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}
    
    tree_1e1a1c={'s1':['d',None, None, 0.1,None,None,None],
      's2':['e',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['a',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s4',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'s2','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}
    
    tree_2a1c={'s1':['c',None, None, 0.1,None,None,None],
      's2':['e',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['a',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s4',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'s2','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}
    
    tree_2c={'s1':['b',None, None, 0.1,None,None,None],
      's2':['e',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s4',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','s1'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'s2','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}
    
    a1=tree_prior(insert_children_in_tree(tree_2e2a))
    a2=tree_prior(tree_3e1a)
    a3=tree_prior(tree_4e)
    a4=tree_prior(tree_2e2c)
    a5=tree_prior(tree_1e1a1c)
    a6=tree_prior(tree_2a1c)
    a7=tree_prior(tree_2c)
    n1=6
    n2=4
    n3=1
    n4=6
    n5=12
    n6=6
    n7=3
    
    print 'tree_2e2a', a1,n1,a1*n1
    print 'tree_3e1a', a2,n2,a2*n2
    print 'tree_4e', a3,n3,a3*n3
    print 'tree_2e2c', a4,n4,a4*n4
    print 'tree_1e1a1c', a5,n5,a5*n5
    print 'tree_2a1c', a6,n6,a6*n6
    print 'tree_2c', a7,n7,a7*n7
    
    print 'total', a1+a2+a3+a4+a5+a6+a7, a1*n1+a2*n2+a3*n3+a4*n4+a5*n5+a6*n6+a7*n7
    
    
    
    