from tree_plotting import plot_as_directed_graph, plot_graph
from Rproposal_admix import addadmix, deladmix
from Rproposal_regraft import make_regraft
from Rproposal_rescale import rescale
from Rtree_operations import create_trivial_tree, make_consistency_checks, get_number_of_admixes, get_trivial_nodes
from numpy.random import choice
from time import sleep as wait
from MCMC import initialize_prior_as_posterior, basic_chain
from summary import s_no_admixes, s_branch_length, s_variable
from meta_proposal import basic_meta_proposal
from generate_prior_trees import generate_admix_topology
from prior import tree_prior
from tree_statistics import unique_identifier
from math import exp
from collections import Counter
from scipy.optimize import brentq

def _get_new_nodes(i,k):
    if k==2:
        return ['x'+str(i)+a for a in ['a','b']]
    else:
        return 'x'+str(i)

def topological_support(start_tree, n=10000, nodes=None):
    tree=start_tree
    score=0
    for i in range(n):
        prop_index=choice(3,1)
        if prop_index==0 and get_number_of_admixes(tree)>0:
            new_tree=deladmix(tree)[0]
            score-=1
        elif prop_index==1:
            new_tree=make_regraft(tree, _get_new_nodes(i, prop_index))[0]
        elif prop_index==2:
            new_tree=addadmix(tree, _get_new_nodes(i, prop_index))[0]
            score+=1
        else:
            new_tree=tree
        consistent, information = make_consistency_checks(new_tree,nodes)
        if consistent:
            if get_number_of_admixes(new_tree)<50:
                tree=new_tree
        else:
            print information
            print 'last_good_tree', tree
            print 'new_bad_tree', new_tree
            plot_as_directed_graph(tree, drawing_name='before.png')
            wait(1)
            plot_as_directed_graph(new_tree, drawing_name='after.png')
            wait(1)
        if i%1000==0:
            plot_as_directed_graph(tree, drawing_name='number_'+str(i)+'K.png')
            wait(1)
    plot_as_directed_graph(tree, drawing_name='final_tree.png')
    wait(1)
    return score

def proposal_support(start_tree, n=10000, nodes=None):
    tree=start_tree
    score=0
    for i in range(n):
        prop_index=choice(4,1)
        if prop_index==0 and get_number_of_admixes(tree)>0:
            new_tree=deladmix(tree)[0]
            score-=1
        elif prop_index==1:
            new_tree=make_regraft(tree, _get_new_nodes(i, prop_index))[0]
        elif prop_index==2:
            new_tree=addadmix(tree, _get_new_nodes(i, prop_index))[0]
            score+=1
        elif prop_index==3:
            new_tree=rescale(tree)[0]
        else:
            new_tree=tree
        consistent, information = make_consistency_checks(new_tree,nodes)
        if consistent:
            if get_number_of_admixes(new_tree)<50:
                tree=new_tree
        else:
            print information
            print 'last_good_tree', tree
            print 'new_bad_tree', new_tree
            plot_as_directed_graph(tree, drawing_name='before.png')
            wait(1)
            plot_as_directed_graph(new_tree, drawing_name='after.png')
            wait(1)
            break
        if i%1000==0:
            print tree
            plot_as_directed_graph(tree, drawing_name='number_'+str(i)+'K.png')
            wait(1)
    plot_as_directed_graph(tree, drawing_name='final_tree.png')
    wait(1)
    return score

def test_prior_model(start_tree, sim_length=100000):
    posterior=initialize_prior_as_posterior()
    summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    proposal=basic_meta_proposal()
    sample_verbose_scheme={'posterior':(1,1),
                           'branch_length':(10,1),
                           'mhr':(1,1),
                           'no_admixes':(1,1)}
    basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=1, i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=True)
    
def test_topological_prior_density(n,k,sim_length):
    dictionary_of_probabilities={}
    list_of_simulated_trees=[]
    for _ in xrange(sim_length):
        tree=generate_admix_topology(n,k)
        unique_id=unique_identifier(tree)
        list_of_simulated_trees.append(unique_id)
        if unique_id not in dictionary_of_probabilities:
            dictionary_of_probabilities[unique_id]=exp(tree_prior(tree))
    ad,ad2=dictionary_of_probabilities, Counter(list_of_simulated_trees)
    for key, val in ad.items():
        print val, 'vs', float(ad2[key])/sim_length,': ', key
    print 'Total number of simulations', sim_length
    print 'Unique sampled trees', len(ad2)
    print 'Simulated total probability', sum(ad.values())
    n=float(sim_length)
    c=float(len(ad2))
    def f(m):
        return m*(1-((m-1)/m)**n)-c
    maxval=c
    while f(maxval)<0:
        maxval*=2
    ex=brentq(f, c, maxval)
    print 'Expected number of unique', ex
    print 'Expected total probability', c/ex
    
    

if __name__=='__main__':
    s_tree=create_trivial_tree(15)
    print #proposal_support(s_tree, nodes= get_trivial_nodes(15))
    #plot_as_directed_graph(s_tree)
    #wait(1)
    test_topological_prior_density(4,2, 10000)

    
    