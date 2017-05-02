from Rproposal_admix import addadmix, deladmix
from Rproposal_regraft import make_regraft
from Rproposal_rescale import rescale
from Rtree_operations import (create_trivial_tree, make_consistency_checks, get_number_of_admixes, get_trivial_nodes, convert_to_vector,pretty_print, pretty_string,
                              get_number_of_ghost_populations, get_max_distance_to_root, get_min_distance_to_root, get_average_distance_to_root, get_no_leaves)
from numpy.random import choice
from time import sleep as wait
from MCMC import basic_chain
from MCMCMC import MCMCMC
from temperature_scheme import fixed_geometrical
from posterior import initialize_prior_as_posterior, initialize_trivial_posterior, initialize_posterior
from summary import s_no_admixes, s_total_branch_length, s_variable, s_posterior, s_average_branch_length, s_total_branch_length, s_basic_tree_statistics
from meta_proposal import basic_meta_proposal, no_admix_proposal, adaptive_proposal, adaptive_proposal_no_admix
from generate_prior_trees import generate_admix_topology, generate_phylogeny
from prior import prior, topological_prior
from tree_statistics import unique_identifier
from math import exp
from collections import Counter
from scipy.optimize import brentq
from analyse_results import save_to_csv
from csv import writer
from trivial_mcmc import Trivial_Summary, trivial_proposal
from Rtree_to_covariance_matrix import make_covariance
from pathos.multiprocessing import Pool
from scipy.stats import geom, wishart
from likelihood import n_mark
import os



def _get_new_nodes(i,k):
    if k==2:
        return ['x'+str(i)+a for a in ['a','b']]
    else:
        return 'x'+str(i)

def topological_support(start_tree, n=10000, nodes=None):
    from tree_plotting import plot_as_directed_graph, plot_graph, pretty_print, pretty_string
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
    proposal=adaptive_proposal() #basic_meta_proposal()
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme['posterior']=(1,1000)
    final_tree,final_posterior, results,_= basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=int(thinning_coef+sim_length/60000), i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=True)
    print results
    save_to_csv(results, summaries)
    return results

def test_prior_model_several_chains(start_trees, sim_length=100000, summaries=None, thinning_coef=1):
    posterior=initialize_prior_as_posterior()
    if summaries is None:
        summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    proposal=basic_meta_proposal()
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    p=Pool(len(start_trees))
    def func(nstart_tree):
        n,start_tree=nstart_tree
        final_tree,final_posterior, results,_= basic_chain(start_tree, summaries, posterior, 
                    proposal, post=None, N=sim_length, 
                    sample_verbose_scheme=sample_verbose_scheme, 
                    overall_thinning=int(thinning_coef+sim_length/60000), i_start_from=0, 
                    temperature=1.0, proposal_update=None,
                    check_trees=True)
        save_to_csv(results, summaries, filename='results_'+str(n+1)+'csv', origin_layer=(n+1,1))
    p.map(func, enumerate(start_trees))

def test_prior_model_no_admixes(start_tree, sim_length=100000, summaries=None, thinning_coef=1):
    posterior=initialize_prior_as_posterior()
    if summaries is None:
        summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    proposal=adaptive_proposal_no_admix()
    #proposal.props=proposal.props[2:] #a little hack under the hood
    #proposal.params=proposal.params[2:] #a little hack under the hood.
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    final_tree,final_posterior, results,_= basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=int(thinning_coef+sim_length/60000), i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=False)
    save_to_csv(results, summaries)
    return results

def test_posterior_model(true_tree=None, start_tree=None, sim_length=100000, summaries=None, thinning_coef=19, admixtures_of_true_tree=None, no_leaves_true_tree=4, filename='results.csv', sim_from_wishart=False, wishart_df=None):
    if true_tree is None:
        if admixtures_of_true_tree is None:
            admixtures_of_true_tree=geom.rvs(p=0.5)-1
        true_tree=generate_phylogeny(no_leaves_true_tree, admixtures_of_true_tree)
    else:
        no_leaves_true_tree=get_no_leaves(true_tree)
        admixtures_of_true_tree=get_number_of_admixes(true_tree)
        
    m=make_covariance(true_tree, get_trivial_nodes(no_leaves_true_tree))
    if start_tree is None:
        start_tree=true_tree
    if wishart_df is None:
        wishart_df=n_mark(m)
    if sim_from_wishart:
        r=m.shape[0]
        print m
        m=wishart.rvs(df=r*wishart_df-1, scale=m/(r*wishart_df))
        print m
    posterior=initialize_posterior(m,wishart_df)
    if summaries is None:
        summaries=[s_posterior(), s_variable('mhr'), s_no_admixes()]
    proposal=adaptive_proposal()
    #proposal.props=proposal.props[2:] #a little hack under the hood
    #proposal.params=proposal.params[2:] #a little hack under the hood.
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme['posterior']=(1,1)
    final_tree,final_posterior, results,_= basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=int(max(thinning_coef,sim_length/60000)), i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=False)
    save_to_csv(results, summaries, filename= filename)
    return true_tree

def test_posterior_model_multichain(true_tree=None, start_tree=None, sim_lengths=[250]*800, summaries=None, thinning_coef=1, admixtures_of_true_tree=None, no_leaves_true_tree=4, wishart_df=None, sim_from_wishart=False, no_chains=4):
    if true_tree is None:
        if admixtures_of_true_tree is None:
            admixtures_of_true_tree=geom.rvs(p=0.5)-1
        true_tree=generate_phylogeny(no_leaves_true_tree, admixtures_of_true_tree)
    else:
        no_leaves_true_tree=get_no_leaves(true_tree)
        admixtures_of_true_tree=get_number_of_admixes(true_tree)
        
    m=make_covariance(true_tree, get_trivial_nodes(no_leaves_true_tree))
    if start_tree is None:
        start_tree=true_tree
    if wishart_df is None:
        wishart_df=n_mark(m)
    if sim_from_wishart:
        r=m.shape[0]
        print m
        m=wishart.rvs(df=r*wishart_df-1, scale=m/(r*wishart_df))
        print m
    posterior=initialize_posterior(m,wishart_df)
    if summaries is None:
        summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    proposal=basic_meta_proposal()
    #proposal.props=proposal.props[2:] #a little hack under the hood
    #proposal.params=proposal.params[2:] #a little hack under the hood.
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme_first=deepcopy(sample_verbose_scheme)
    sample_verbose_scheme_first['posterior']=(1,1)
    final_tree,final_posterior, results,_=MCMCMC(starting_trees=[deepcopy(start_tree) for _ in range(no_chains)], 
               posterior_function= posterior,
               summaries=summaries, 
               temperature_scheme=fixed_geometrical(10.0,no_chains), 
               printing_schemes=[sample_verbose_scheme_first]+[sample_verbose_scheme for _ in range(no_chains-1)], 
               iteration_scheme=sim_lengths, 
               overall_thinnings=int(thinning_coef+sim_length/60000), 
               proposal_scheme= [proposal_function for _ in range(n)], 
               cores=no_chains, 
               no_chains=no_chains)
    
    basic_chain(start_tree, summaries, posterior, 
                proposal, post=None, N=sim_length, 
                sample_verbose_scheme=sample_verbose_scheme, 
                overall_thinning=int(thinning_coef+sim_length/60000), i_start_from=0, 
                temperature=1.0, proposal_update=None,
                check_trees=False)
    save_to_csv(results, summaries)
    return true_tree
    
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
        tup=(pks['source_key'],pks['source_branch'],pks['sink_key'],pks['sink_branch'], pks['sink_new_branch'])
        if tup in dic:
            assert dic[tup][0]==unique_id, 'There is not a unique mapping from sink source to unique id'
            assert abs(1.0/dic[tup][1]-pks['forward_choices'])<1e-8, 'There was not a unqiue number of forward choices for '+str(tup)+' : '+str(pks['forward_choices'])+'><' +str(1.0/dic[tup][1])
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
        tup=(pks['remove_key'],pks['remove_branch'])#,pks['sink_key'],pks['sink_branch'])
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

def marginalize_out_data_in_posterior(no_leaves, no_trees=100, nsim=50000, wishart_df=None, prefix='', dest_folder=''):
    summaries=[s_posterior(), 
           s_no_admixes(),
           s_average_branch_length(),
           s_total_branch_length(),
           s_basic_tree_statistics(get_number_of_ghost_populations, 'ghost_pops', output='integer'),
           s_basic_tree_statistics(get_max_distance_to_root, 'max_root'),
           s_basic_tree_statistics(get_min_distance_to_root, 'min_root'),
           s_basic_tree_statistics(get_average_distance_to_root, 'average_root')]+[s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    
    for i in xrange(no_trees):
        result_file=os.path.join(dest_folder,'results_'+prefix+str(i+1)+'.csv')
        test_posterior_model(thinning_coef=49, summaries= summaries, filename= result_file, sim_from_wishart=True,sim_length=nsim, wishart_df=wishart_df, no_leaves_true_tree=no_leaves)
    
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
    #s_tree=create_trivial_tree(4)
    
    marginalize_out_data_in_posterior(4, no_trees=50, nsim=100000, wishart_df=10)
    #marginalize_out_data_in_posterior(4, no_trees=50, nsim=100000, wishart_df=100)
    #marginalize_out_data_in_posterior(4, no_trees=50, nsim=100000, wishart_df=1000)
    #topological_support(s_tree)
    #print test_prior_model(s_tree, 100000)
     #proposal_support(s_tree, nodes= get_trivial_nodes(15))
    #plot_as_directed_graph(s_tree)
    #wait(1)
    #print test_topological_prior_density(3,3, 500000)
    #from Rcatalogue_of_trees import tree_good, tree_one_admixture, tree_minimal
    
    #tree2=addadmix(addadmix(addadmix(addadmix(tree_good)[0])[0])[0])[0]
    #print proper_proposals(tree_good,100000)[1]
    #print check_destinations(tree_minimal, 20000)
    #print check_predestinations(tree2, 100000)

    
    