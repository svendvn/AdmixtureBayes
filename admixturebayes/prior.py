from scipy.stats import expon, geom, pareto, nbinom
from Rtree_operations import (get_all_branch_lengths, get_all_admixture_proportions, get_number_of_admixes, get_number_of_leaves, 
get_leaf_keys,get_destination_of_lineages, get_categories, get_parent_of_branch, propagate_married, propagate_admixtures)
from math import log, factorial,exp
from scipy.special import binom 
import linear_distribution
from uniform_topological_prior import uniform_prior, uniform_topological_prior_function
from Rtree_to_covariance_matrix import get_admixtured_populations

def calculate_branch_prior(branches, n):
    #if all((b<10 for b in branches)):
    #    return -len(branches)*log(10)
    #else:
    #    return -float('inf')
    #return 0.0
    #return sum((pareto.logpdf(br,2, scale=0.01) for br in branches))
    rate=float(2*n-2)/len(branches)
    n2k=len(branches)
    d=float(n2k)/float(n)
    
    return -sum(branches)*rate+log(rate)*len(branches)
    return -sum(branches)
    #return sum(map(expon.logpdf, branches))

def illegal_admixtures(unadmixed_populations, tree):
    admixed_populations= get_admixtured_populations(tree)
    if set(admixed_populations).intersection(unadmixed_populations):
        return True
    return False

def calculate_add_prior(add, rate=1):
    '''The rate is the mean.'''
    return -add/float(rate)


def prior(x, p=0.5, use_skewed_distr=False, pks={}, use_uniform_prior=False, unadmixed_populations=[], r=0, add_prior_rate=1):
    tree, add=x
    no_leaves=get_number_of_leaves(tree)
    admixtures=get_all_admixture_proportions(tree)
    if not all(prop>=0 and prop <=1 for prop in admixtures):
        return -float('inf')
    branches=get_all_branch_lengths(tree)
    if not all(branch>=0 for branch in branches):
        return -float('inf')
    branch_prior=calculate_branch_prior(branches, no_leaves)
    no_admix_prior=no_admixes(p, len(admixtures), r=r)
    if use_skewed_distr:
        admix_prop_prior=linear_admixture_proportions(admixtures)
    else:
        admix_prop_prior=0
    if use_uniform_prior:
        top_prior=uniform_topological_prior_function(tree)
    else:
        top_prior=topological_prior(tree)
    if unadmixed_populations:
        if illegal_admixtures(unadmixed_populations, tree):
            return -float('inf')
    add_prior=calculate_add_prior(add, add_prior_rate)
    logsum=branch_prior+no_admix_prior+admix_prop_prior+top_prior+add_prior
    pks['branch_prior']= branch_prior
    pks['no_admix_prior']=no_admix_prior
    pks['admix_prop_prior']=admix_prop_prior
    pks['top_prior']= top_prior
    return logsum

def linear_admixture_proportions(admixtures):
    return sum((linear_distribution.logpdf(admixture) for admixture in admixtures))

def no_admixes(p, admixes, hard_cutoff=20, r=0):
    if admixes>hard_cutoff:
        return -float('inf')
    if r>1:
        if hard_cutoff is None:
            return nbinom.logpmf(admixes,n=r, p=1.0-p)
        else:
            return nbinom.logpmf(admixes, n=r, p=1.0 - p) - nbinom.logcdf(hard_cutoff ,n=r, p= 1.0 - p)
    else:
        if hard_cutoff is None:
            return geom.logpmf(admixes+1, 1.0-p)

        return geom.logpmf(admixes+1, 1.0-p)-geom.logcdf(hard_cutoff+1, 1.0-p)
    
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
    '''
    calculates the probability of the selection of the number of totally free, halfly free and admixture events given the availabe number of each type.
    '''
    #print total_admixtures, total_occupied, total_unoccupied
    #print sampled_admixtures, sampled_occupied, sampled_unoccupied
    
    logmhypergeom= (log(binom(total_unoccupied,sampled_unoccupied))+
            log(binom(total_occupied,sampled_occupied))+
            log(binom(total_admixtures, sampled_admixtures))-
            log(binom(total_unoccupied+total_occupied+total_admixtures,
                      sampled_unoccupied+sampled_occupied+sampled_admixtures)))
    
    return logmhypergeom
    
def _totally_free_selection(total_pairs,sampled_pairs, sampled_ones):
    '''
    calculates the probability that the distribution of 1's and 2's, and 0's are as specified in sampled_pairs, sampled_ones. 
    A 1 means that one of the completely free coalescence only one coalescence holster has been taken
    '''
    denom=log(binom(total_pairs*2,sampled_pairs*2+sampled_ones))
    num=log(2)*sampled_ones
    num+=log(factorial(total_pairs))
    num-=sum(log(factorial(x)) for x in [total_pairs-sampled_ones-sampled_pairs,
                                         sampled_pairs,
                                         sampled_ones])
    return num-denom

def _awaited_selection(positions, taken):
    '''
    Calculates the probability that exactly the of the number of available waiting coalescence nodes (=positions), the exact subset of size 'taken' is taken out from them. 
    '''
    return -log(binom(positions,taken))

def _assigned_selection(total, pairs, unlabeleld_singles, labelled_singles, admixtures):
    '''
    Given the number of each type of parent key, this function calculates the probability of the way it has been laid out.
    '''
    return -(log(factorial(total))-
             log(factorial(pairs))-
             log(factorial(unlabeleld_singles))-
             log(factorial(admixtures))-
             log(2)*pairs)
    
def _get_number_of_admixture_branches(lineages):
    return sum(1 for key, branch in lineages if branch==1)

def _double_band_selection(no_danger_pairs, M,N):
    '''
    calculated that none of the pairs no_danger_pairs makes a double band when there are a total of N destinations and M of these are totally free coalescence nodes
    '''
    if N<M:
        return 0
    if no_danger_pairs==0:
        return 1
    if M<2:
        return 1
    res=0
    p1=float(M*(M-2))/N/(N-1)
    next_step=_double_band_selection(no_danger_pairs-1, M-4, N-2)
    res+=p1*next_step
    if N>2:
        p2=float(M*(N-M))/N/(N-1)
        p3=float((N-M)*M)/N/(N-1)
        next_step2=_double_band_selection(no_danger_pairs-1, M-2, N-2)
        res+=(p2+p3)*next_step2
        if N>3:
            p4=float((N-M)*(N-M-1))/N/(N-1)
            next_step3=_double_band_selection(no_danger_pairs-1, M, N-2)
            res+=next_step3*p4
    
    return res

def topological_prior(tree):
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
    
    
    #return sames, single_coalescences, admixtures, waiting_coalescences,awaited_coalescences,still_on_hold,p1,p2,p3,p4, exp(p1+p2+p3+p4)
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
    p_other=0
    if no_free_admixtures==1 and no_totally_free_coalescences==1 and no_coalescences_on_hold==1 and total_chosen==1:#in this case we have to choose the admixture, and there are no choice left
        return 0
    if no_free_admixtures==1 and no_totally_free_coalescences==2 and no_coalescences_on_hold==0 and total_chosen==2:#here the choices are either 1admixture+1coalescence or 1coalescence+1admixture and there are only two options, 
        return -log(2)
    p_no_double_bands=_double_band_selection(no_admixture_pairs, 
                                                 2*no_totally_free_coalescences, 
                                                 2*no_totally_free_coalescences+no_coalescences_on_hold+no_free_admixtures)
    logp_free_coalescences=0
    if no_totally_free_coalescences>=total_chosen:
        logp_free_coalescences+=_events_selection(sampled_unoccupied=total_chosen, 
                            sampled_occupied=0, 
                            sampled_admixtures=0,
                            total_unoccupied=2*no_totally_free_coalescences,
                            total_occupied=no_coalescences_on_hold,
                            total_admixtures=no_free_admixtures)
        logp_free_coalescences+=_totally_free_selection(no_totally_free_coalescences, 0, total_chosen)
        p_admix_minus=log(1.0-exp(logp_free_coalescences)-(1.0-p_no_double_bands))
    else:
        p_admix_minus=log(p_no_double_bands)
    #print exp(logp_free_coalescences), p_no_double_bands
    #print total_chosen, no_totally_free_coalescences, no_free_admixtures, no_coalescences_on_hold
    
    #if no_admixture_pairs>0:
    #val=_double_band_selection(no_admixture_pairs, 
    #                                             2*no_totally_free_coalescences, 
    #                                             2*no_totally_free_coalescences+no_coalescences_on_hold+no_free_admixtures)
#     if val<=0:
#         print val
#         val=float('inf')
#         print no_admixture_pairs, 2*no_totally_free_coalescences, 2*no_totally_free_coalescences+no_coalescences_on_hold+no_free_admixtures

    p3=_awaited_selection(no_coalescences_on_hold,no_awaited_coalescences)
    p4=_assigned_selection(no_ready_lineages, no_sames, no_waiting_coalescences, no_awaited_coalescences, no_admixtures)
    return p1+p2+p3+p4-p_admix_minus

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
    print exp(get_admixture_extra_tops(3, 2))
    from Rcatalogue_of_trees import tree_good, tree_clean, tree_one_admixture, tree_two_admixture
    from Rtree_operations import insert_children_in_tree
    
    print _double_band_selection(2,4,4)
    

    
    #print prior(tree_good)
    
    #print prior(insert_children_in_tree(tree_one_admixture))
    #print prior(insert_children_in_tree(tree_two_admixture))
    
    #res=tree_prior(insert_children_in_tree(tree_clean))
    #print res, exp(res)

    
    from generate_prior_trees import generate_admix_topology
    from tree_plotting import pretty_print
    
    print 'ttt'
    for i in range(25):
        print exp(no_admixes(0.5,i,20))
    
    tree_trouble={'s3': ['n1', None, None, 0.13, None, None, None], 's2': ['n2', None, None, 0.14, None, None, None], 's1': ['n1', None, None, 0.14, None, None, None], 'a3': ['r', 'r', 0.6, 0.14, 0.13, 'n2', None], 'n1': ['n2', None, None, 0.13, None, 's1', 's3'], 'n2': ['a3', None, None, 0.12, None, 's2', 'n1']}
    #print tree_prior(tree_trouble)
    
    for i in range(200):
        tree=generate_admix_topology(3,1)
        #print tree
        #pretty_print(tree)
        #print exp(prior((tree,0))),
        
    import sys
    sys.exit()
        
    
    
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
    
    
    
    