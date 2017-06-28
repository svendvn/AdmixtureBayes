import os 
#os.environ["OPENBLAS_NUM_THREADS"] = "1"
#os.environ["MKL_NUM_THREADS"] = "1"
from MCMC import basic_chain
from pathos.multiprocessing import Pool, freeze_support
import pandas as pd
from multiprocessing_helpers import basic_chain_pool
from numpy.random import choice, random
from math import exp
from itertools import chain


def _basic_chain_unpacker(args):
    return basic_chain(*args)

def MCMCMC(starting_trees, 
           posterior_function,
           summaries,
           temperature_scheme, 
           printing_schemes, 
           iteration_scheme, 
           overall_thinnings, 
           proposal_scheme, 
           cores=4,
           no_chains=None,
           numpy_seeds=None,
           multiplier= None):
    '''
    this function runs a MC3 using the basic_chain_unpacker. Let no_chains=number of chains. The inputs are
        starting_trees: a list of one or more trees that the chains should be started with
        posterior_function: one 'initialized'(see MCMC.py) unnormalized posterior function, that should be simulated from.
        summaries: a list of instances of realization of concrete classes of the superclass Summary. It is closely linked to printing_scheme, 
                   because they are only meaningful if specified in that.
        temperature_scheme: one instance of a class that has the functions:
                                    - get_temp(i): returns the the temperature for the i'th chain
                                    - update_temp(permut): updates the temperatures for each chain using the permutation, permut
        printing_scheme: a list of either one or no_chains dictionaries, where each dictionary gives the sample_verbose_scheme of basic_chain
        iteration_scheme: a list that sums to the total number of iterations, where each entry is the number of MH-steps before the next flipping
        overall_thinnings: an integer indicating how much should be skipped before each summary statistics is calculated
        proposal_scheme: a list of instances of classes that handles proposals and updates of said proposals. It has the basic functions
                                    - prop(x,pks): proposes the next tree(and returns proposal densities and statistics in pks)
                                    - adapt(mhr): updates the proposals based on the mhr ratio, mhr
                                    - extract_new_values(): after two chains have changed position in the flipping, this function gets the information that should be changed_later
                                    - wear_new_values(information): the new values from extract_new_values replaces the old values.
        #ps:     - a list of doubles in the interval [0,1] indicating the parameter of the geometrically distributed prior on the number of admixture events for each chain. Or:
        #        - a double in the interval [0,1] indicating the same parameter for all geometrically distributed prior on the number of admixture events.
                                    
    '''
    
    if no_chains is None:
        no_chains=cores

    #if isinstance(ps, int):
     #   ps=[ps]*no_chains

        
    if len(printing_schemes)==1:
        printing_schemes=[printing_schemes[0]]*no_chains
        

    
    df_result=None
    total_permutation=range(no_chains)
    permuts=[]
    
    if numpy_seeds is None:
        numpy_seeds=[None]*no_chains
    pool = basic_chain_pool(summaries, posterior_function, proposal_scheme, numpy_seeds)
        
    
    xs = starting_trees
    posteriors = [posterior_function(x) for x in xs]
    proposal_updates=[proposal.get_exportable_state() for proposal in proposal_scheme]
    
    cum_iterations=0
    for no_iterations in iteration_scheme:
        #letting each chain run for no_iterations:
        iteration_object=_pack_everything(xs, posteriors, temperature_scheme, printing_schemes, overall_thinnings, no_iterations, cum_iterations, proposal_updates, multiplier)
        if no_chains==1:#debugging purposes
            new_state=[_basic_chain_unpacker(iteration_object.next())]
        else:
            new_state = pool.order_calculation(iteration_object)
        xs, posteriors, df_add,proposal_updates = _unpack_everything(new_state, summaries, total_permutation)
        df_result=_update_results(df_result, df_add)
        
        #making the mc3 flips and updating:
        xs, posteriors, permut, proposal_updates=flipping(xs, posteriors, temperature_scheme, proposal_updates) #trees, posteriors, range(len(trees)),[None]*len(trees)#
        permuts.append(permut)
        temperature_scheme.update_temps(permut)
        #proposal_updates=_handle_flipping_of_proposals(proposal_updates, permut)
        total_permutation=_update_permutation(total_permutation, permut)
        cum_iterations+=no_iterations
        
    pool.terminate()
        
    return df_result, permuts
        
def _update_permutation(config, permut):
    return [config[n] for n in permut]
    
def flipping(xs, posteriors, temperature_scheme, proposal_updates):
    n=len(xs)
    step_permutation=range(n)
    count=0
    for _ in xrange(40):
        i,j = choice(n,2,False)
        post_i,post_j=posteriors[i],posteriors[j]
        temp_i,temp_j=temperature_scheme.get_temp(i), temperature_scheme.get_temp(j)
        logalpha=-(post_i[0]+post_i[1]-post_j[0]-post_j[1])*(1.0/temp_i-1.0/temp_j)
        #print alpha
        if logalpha>0 or random() < exp(logalpha):
            count+=1
            #if i==0 or j==0:
                #print 'FLIP!', count
                #print temp_i, post_i
                #print temp_j, post_j
            step_permutation[i],step_permutation[j]=step_permutation[j],step_permutation[i]
            posteriors[j],posteriors[i]=post_i,post_j
            xs[i],xs[j]=xs[j],xs[i]
            proposal_updates[i],proposal_updates[j]=proposal_updates[j],proposal_updates[i]
    #print step_permutation
    return xs, posteriors, step_permutation, proposal_updates

def _update_results(df_result, df_add):
    if df_result is None:
        df_result = df_add
    else:
        df_result = pd.concat([df_result, df_add])
    return df_result

##should match basic_chain_class.run_chain(start_tree, post, N, sample_verbose_scheme, overall_thinning, i_start_from, temperature, proposal_update)
def _pack_everything(xs, posteriors, temperature_scheme,printing_schemes,overall_thinnings,no_iterations,cum_iterations, proposal_updates=None, multiplier=None):
    return ([x, 
             posterior,
             no_iterations,
             printing_scheme,
             overall_thinnings,
             cum_iterations,
             temperature_scheme.get_temp(i),
             proposal_update,
             multiplier] for i,(x,posterior,printing_scheme,proposal_update) in enumerate(zip(xs,posteriors,printing_schemes,proposal_updates)))

def _unpack_everything(new_state, summaries, total_permutation):
    xs,posteriors, summs, proposal_updates = zip(*new_state)
    list_of_smaller_data_frames=[]
    for summ_data, n, i in zip(summs, total_permutation, range(len(total_permutation))):
        iter_chain=chain((('iteration', summ_data[0]),),
                         ((summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,summ_data[1:])))
        df=pd.DataFrame.from_items(iter_chain)
        df['origin']=n
        df['layer']=i
        list_of_smaller_data_frames.append(df)
    df=pd.concat(list_of_smaller_data_frames)
    return list(xs), list(posteriors), df, list(proposal_updates)

def _handle_flipping_of_proposals(proposal_scheme, permut):
    if permut==range(len(permut)):
        return proposal_scheme
    places_that_change_position={n:m for n,m in enumerate(permut) if n!=m}
    emigrating_values=[]
    for position in places_that_change_position.values():
        emigrating_values.append(proposal_scheme[position].extract_new_values())
    for position in places_that_change_position: ##HUSK AT HOLDE STYR PAA PERMUTATIONERNE!
        proposal_scheme[position].wear_new_values(emigrating_values.pop(0))
    return proposal_scheme

def run_test():
    from Rtree_operations import get_trivial_nodes, create_trivial_tree,get_number_of_ghost_populations,get_max_distance_to_root,get_min_distance_to_root,get_average_distance_to_root
    from posterior import initialize_prior_as_posterior, initialize_posterior
    from meta_proposal import basic_meta_proposal
    from copy import deepcopy
    from Rtree_to_covariance_matrix import make_covariance
    
    
    
    N=3
    true_tree=create_trivial_tree(N)
    proposal_function= basic_meta_proposal()
    post_fun=initialize_posterior(make_covariance(true_tree))
    tree= create_trivial_tree(N)
    
    n=6
    import summary
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(get_average_distance_to_root, 'average_root'),
               summary.s_variable('proposal_type', output='string')]
    
    from temperature_scheme import fixed_geometrical
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme['posterior']=(1,100)
    #sample_verbose_scheme['min_root']=(1,100)
    
    ad=MCMCMC(starting_trees=[deepcopy(tree) for _ in range(n)], 
               posterior_function= post_fun,
               summaries=summaries, 
               temperature_scheme=fixed_geometrical(10.0,n), 
               printing_schemes=[sample_verbose_scheme for _ in range(n)], 
               iteration_scheme=[40]*200, 
               overall_thinnings=5, 
               proposal_scheme= [proposal_function for _ in range(n)], 
               cores=n, 
               no_chains=n)
    
    ad[0].to_csv(path_or_buf='findme.csv')
    print set(map(tuple,ad[1]))
    return ad

if __name__=='__main__':
        
    run_test()
    
#     from scipy.stats import poisson
#     from likelihood import likelihood
#     from numpy.random import random
#     from math import exp
#     from tree_operations import get_number_of_admixes
#     from summary import *
#     import itertools
#     
#     from Proposal_Function import prop_flat
#         
#         
#     from tree_operations import make_flat_list_no_admix
#         
#     from numpy import diag
#     N=10
#     tree_flatter_list=make_flat_list_no_admix(N)
#     nodes=["s"+str(i) for i in range(1,N+1)]
#     emp_cov=diag([0.5]*N)
#     emp_cov[2,1]=emp_cov[1,2]=0.2
#     
#     x=tree_flatter_list
#     posterior_function=initialize_posterior(emp_cov)
#     summaries=[s_variable('posterior'), s_variable('mhr'), s_branch_length()]
#     sample_verbose_scheme={'posterior':(1,100),
#                            'branch_length':(5,0),
#                            'mhr':(1,0)}
#     n=2
#     def run():
#         ad=MCMCMC(starting_trees=[x for _ in range(n)], 
#                posterior_function= posterior_function,
#                summaries=summaries, 
#                temperature_scheme=fixed_geometrical(10.0,n), 
#                printing_schemes=[sample_verbose_scheme for _ in range(n)], 
#                iteration_scheme=[10]*100, 
#                overall_thinnings=5, 
#                proposal_scheme= [prop_flat for _ in range(n)], 
#                cores=n, 
#                no_chains=n)
#     #run()
#     import cProfile
#     cProfile.run('run()')
    
