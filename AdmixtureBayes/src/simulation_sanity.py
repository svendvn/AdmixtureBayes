from tree_plotting import plot_as_directed_graph, plot_graph, pretty_print, pretty_string
from Rproposal_admix import addadmix, deladmix
from Rproposal_regraft import make_regraft
from Rproposal_rescale import rescale
from Rtree_operations import create_trivial_tree, make_consistency_checks, get_number_of_admixes, get_trivial_nodes, convert_to_vector
from numpy.random import choice
from time import sleep as wait
from MCMC import basic_chain
from posterior import initialize_prior_as_posterior, initialize_trivial_posterior
from summary import s_no_admixes, s_total_branch_length, s_variable
from meta_proposal import basic_meta_proposal
from generate_prior_trees import generate_admix_topology
from prior import prior, topological_prior
from tree_statistics import unique_identifier
from math import exp
from collections import Counter
from scipy.optimize import brentq
from analyse_results import save_to_csv
from csv import writer
from trivial_mcmc import Trivial_Summary, trivial_proposal


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

def test_prior_model(start_tree, sim_length=100000, summaries=None, thinning_coef=1):
    posterior=initialize_prior_as_posterior()
    if summaries is None:
        summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    proposal=basic_meta_proposal()
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    final_tree,final_posterior, results,_= basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=int(thinning_coef+sim_length/60000), i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=True)
    print results
    save_to_csv(results, summaries)
    return results

def test_prior_model_no_admixes(start_tree, sim_length=100000, summaries=None):
    posterior=initialize_prior_as_posterior()
    if summaries is None:
        summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    proposal=basic_meta_proposal()
    proposal.props=proposal.props[2:] #a little hack under the hood
    proposal.params=proposal.params[2:] #a little hack under the hood.
    sample_verbose_scheme={summary.name:(1,1) for summary in summaries}
    final_tree,final_posterior, results,_= basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=1, i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=True)
    save_to_csv(results, summaries)
    return results
    
def test_topological_prior_density(n,k,sim_length):
    dictionary_of_probabilities={}
    list_of_simulated_trees=[]
    for _ in xrange(sim_length):
        tree=generate_admix_topology(n,k)
        unique_id=unique_identifier(tree)
        list_of_simulated_trees.append(unique_id)
        if unique_id not in dictionary_of_probabilities:
            print unique_id
            dictionary_of_probabilities[unique_id]=exp(topological_prior(tree))
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
    
def proper_proposals(start_tree, reps=10000):
    points=[]
    densities=[]
    keys=['s1','s2','s3','s4','a','b','c','d','e','f','x','y']
    for _ in xrange(reps):
        pks={}
        proposed,forward, backward = addadmix(start_tree,new_node_names=['x','y'], pks=pks, fixed_sink_source = ('a',0,'s1',0))
        numbers=convert_to_vector(proposed,keys=keys)
        density=pks['forward_density']
        points.append(numbers+[density])

    with open("results2.csv", "wb") as f:
        wr = writer(f)
        wr.writerows(points)
    return 'done'

def check_destinations(tree, reps=10000):
    dic={}
    
    for _ in xrange(reps):
        pks={}
        ntree, _, _ = addadmix(tree, new_node_names=['x','y'],pks=pks)
        unique_id=unique_identifier(ntree)
        tup=(pks['source_key'],pks['source_branch'],pks['sink_key'],pks['sink_branch'])
        if tup in dic:
            assert dic[tup][0]==unique_id, 'There is not a unique mapping from sink source to unique id'
            assert 1.0/dic[tup][1]==pks['forward_choices'], 'There was not a unqiue number of forward choices for '#+str(tup)+' : '+str(pks['forward_choices'])+'><'
            dic[tup]=(unique_id, 1.0/pks['forward_choices'], dic[tup][2]+1)
        else:
            dic[tup]=(unique_id, 1.0/pks['forward_choices'],1)
    res_ret=[]
    total=0
    for key, element in dic.items():
        print key, ': ', element
        total+=element[1]
        res_ret.append(list(key)+list(element))
    print 'TOTAL=', total
    with open("results3.csv", "wb") as f:
        wr = writer(f)
        wr.writerows(res_ret)
    return 'done'
        
def check_predestinations(tree, reps=10000):
    dic={}
    
    for _ in xrange(reps):
        pks={}
        ntree, _, _=deladmix(tree,pks=pks)
        unique_id=unique_identifier(ntree)
        tup=(pks['remove_key'],pks['remove_branch'],pks['sink_key'],pks['sink_branch'])
        if tup in dic:
            assert dic[tup][0]==unique_id, 'There is not a unique mapping from sink source to unique id'
            assert 1.0/dic[tup][1]==pks['forward_choices'], 'There was not a unqiue number of forward choices for '#+str(tup)+' : '+str(pks['forward_choices'])+'><'
            dic[tup]=(unique_id, 1.0/pks['forward_choices'], dic[tup][2]+1)
        else:
            dic[tup]=(unique_id, 1.0/pks['forward_choices'],1)
    res_ret=[]
    total=0
    for key, element in dic.items():
        print key, ': ', element
        total+=element[1]
        res_ret.append(list(key)+list(element))
    print 'TOTAL=', total
    with open("results3.csv", "wb") as f:
        wr = writer(f)
        wr.writerows(res_ret)
    return 'done'     
    
def trivial_simulation(start_val, reps, thinning_coef=1):
    posterior=initialize_trivial_posterior()
    summaries=[Trivial_Summary()]
    proposal=trivial_proposal()
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    final_tree,final_posterior, results,_= basic_chain(start_val, summaries, posterior, 
                proposal, post=None, N=reps, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=int(thinning_coef+reps/60000), i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=False)
    print results
    save_to_csv(results, summaries)
    return results



if __name__=='__main__':
    s_tree=create_trivial_tree(4)
    topological_support(s_tree)
    #print test_prior_model(s_tree, 100000)
     #proposal_support(s_tree, nodes= get_trivial_nodes(15))
    #plot_as_directed_graph(s_tree)
    #wait(1)
    #print test_topological_prior_density(3,3, 500000)
    from Rcatalogue_of_trees import tree_good
    #print proper_proposals(tree_good,100000)[1]
    #print check_predestinations(tree_good, 10000)

    
    