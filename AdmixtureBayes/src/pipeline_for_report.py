import simulation_sanity
import analyse_results
import Rcatalogue_of_trees
import Rtree_operations
import summary
import generate_prior_trees
import trivial_mcmc
import tree_statistics
#import tree_generation_laboratory
import load_data
from numpy import diag, ndarray
import tree_to_data

def run_a():
    n=4
    s_tree=Rtree_operations.create_burled_leaved_tree(n,1)
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
               
    simulation_sanity.test_prior_model(s_tree, 50000, summaries=summaries, thinning_coef=3)
    def max_two(tree):
        if Rtree_operations.get_number_of_admixes(tree)>2:
            return False
        return True
    list_of_summaries=summaries[2:10]
    nsim=100000
    #prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, sim_length=nsim, list_of_summaries=list_of_summaries, skewed_admixture_prior=False)#, thinning_criteria=max_two)
    #analyse_results.save_to_csv([tuple(range(nsim))]+[tuple(prior_distribution[summ.name]) for summ in list_of_summaries], list_of_summaries, filename='sim_prior2.csv', origin_layer=None)
    #analyse_results.generate_summary_csv(summaries)
    #analyse_results.full_analysis(summaries,
    #              trajectories_for_all_temperatures=False,
    #              trajectories_for_all_chains=False,
    #              prior_distribution=prior_distribution)
    
def run_b():
    n=4
    s_tree=Rtree_operations.create_trivial_tree(n)
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]    
    simulation_sanity.test_prior_model_no_admixes(s_tree, 200000, summaries=summaries, thinning_coef=3)
    list_of_summaries=summaries[2:10]
    nsim=100000
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, admixes=0, sim_length=nsim, list_of_summaries=list_of_summaries)
    analyse_results.save_to_csv([tuple(range(nsim))]+[tuple(prior_distribution[summ.name]) for summ in list_of_summaries], list_of_summaries, filename='sim_prior.csv', origin_layer=None)
    analyse_results.generate_summary_csv(summaries)
    #analyse_results.full_analysis(summaries,
    #              trajectories_for_all_temperatures=False,
    #              trajectories_for_all_chains=False,
    #              prior_distribution=prior_distribution)
    
def run_c():
    n=3
    s_trees=[Rtree_operations.create_trivial_tree(n), Rtree_operations.create_burled_leaved_tree(n,1.0), Rtree_operations.create_balanced_tree(n,1.0)]
    summaries=[summary.s_variable('posterior'), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
               
    simulation_sanity.test_prior_model_several_chains(s_trees, 100000, summaries=summaries, thinning_coef=3)
    print 'finished mcmc chains'
    list_of_summaries=summaries[2:10]
    nsim=100000
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, sim_length=nsim, list_of_summaries=list_of_summaries)#, thinning_criteria=max_two)
    analyse_results.save_to_csv([tuple(range(nsim))]+[tuple(prior_distribution[summ.name]) for summ in list_of_summaries], list_of_summaries, filename='sim_prior.csv', origin_layer=None)
    analyse_results.generate_summary_csv(summaries)
    
def run_d(true_tree_as_file=None):
    #true_tree=generate_prior_trees.generate_phylogeny(8,2)
    if true_tree_as_file is None:
        true_tree=tree_statistics.identifier_to_tree_clean('w.w.w.w.w.w.a.a.w-c.w.c.c.w.c.5.0.w.3.2-c.w.w.0.c.4.w-c.w.0.c.3-w.c.1-c.0;0.07-0.974-1.016-0.089-0.81-0.086-1.499-0.052-1.199-2.86-0.403-0.468-0.469-1.348-1.302-1.832-0.288-0.18-0.45-0.922-2.925-3.403;0.388-0.485')
        #true_tree=Rcatalogue_of_trees.tree_good
        s_tree=tree_statistics.identifier_to_tree_clean('w.w.a.w.w.a.a.a.w-c.w.c.c.w.w.c.0.w.w.6.3.2-c.w.w.0.w.c.5.w.w-c.w.0.c.3.w.w-c.w.c.2.0-w.c.1-c.0;0.828-0.21-0.197-0.247-0.568-1.06-0.799-1.162-2.632-2.001-0.45-1.048-0.834-0.469-0.191-2.759-0.871-1.896-0.473-0.019-1.236-0.287-0.179-0.981-0.456-0.91-2.114-3.368;0.655-0.506-0.389-0.23')
        print Rtree_operations.pretty_string(s_tree)
        print Rtree_operations.pretty_string(true_tree)
    else:
        with open(true_tree_as_file, 'r') as f:
            s=f.readline().rstrip()
            true_tree=tree_statistics.identifier_to_tree_clean(s)
            no_leaves= Rtree_operations.get_number_of_leaves(true_tree)
            s_tree= Rtree_operations.create_trivial_tree(no_leaves)
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr', output='double_missing'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_basic_tree_statistics(tree_statistics.get_admixture_proportion_string, 'admixtures', output='string'),
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_basic_tree_statistics(tree_statistics.majority_tree, 'majority_tree', output='string'),
               summary.s_variable('add', output='double'),
               summary.s_variable('sliding_rescale_adap_param', output= 'double_missing'),
               summary.s_variable('cutoff_distance', output= 'double_missing'),
               summary.s_variable('number_of_pieces', output= 'double_missing'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double_missing') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model(true_tree,s_tree, 20000, summaries=summaries, thinning_coef=20, wishart_df= 10000, resimulate_regrafted_branch_length=False)#, 
                                             #admixtures_of_true_tree=2, no_leaves_true_tree=8, rescale_empirical_cov=True)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(r)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree)

def max_dist(x,y):
    signed_diffs= [(abs(xi-yi), xi-yi) for xi,yi in zip(x,y)]
    return sorted(signed_diffs)[-1][1]

def run_analysis_of_proposals():
    #true_tree=generate_prior_trees.generate_phylogeny(8,2)
    true_tree=tree_statistics.identifier_to_tree_clean('w.w.c.w.w.w.2.w-w.w.a.w.w.w.w-w.c.1.w.c.w.w.4-w.c.1.w.w.w-w.c.1.w.w-c.0.w.w-c.w.0-a.w-c.0.w-c.0;0.091-1.665-0.263-0.821-0.058-0.501-0.141-0.868-5.064-0.153-0.372-3.715-1.234-0.913-2.186-0.168-0.542-0.056-2.558-0.324;0.367-0.451')
    true_tree=Rcatalogue_of_trees.tree_good
    s_tree=Rtree_operations.create_trivial_tree(4)
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_basic_tree_statistics(tree_statistics.majority_tree, 'majority_tree', output='string'),
               summary.s_bposterior_difference(lambda x: x[0], 'likelihood_difference'),
               summary.s_bposterior_difference(lambda x: x[1], 'prior_difference'),
               summary.s_bposterior_difference(lambda x: x[2][0], 'branch_prior_difference'),
               summary.s_bposterior_difference(lambda x: x[2][1], 'no_admix_prior_difference'),
               summary.s_bposterior_difference(lambda x: x[2][2], 'adix_prop_prior_difference'),
               summary.s_bposterior_difference(lambda x: x[2][3], 'top_prior_difference'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double_missing') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model(true_tree,true_tree, 20000, summaries=summaries, thinning_coef=2, wishart_df= 1000, resimulate_regrafted_branch_length=False, 
                                             admixtures_of_true_tree=2, no_leaves_true_tree=4, big_posterior=True, rescale_empirical_cov=True)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(r)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree)

    
def run_e(df, out, sap_sim, sap_ana):
    simulation_sanity.marginalize_out_data_in_posterior(4, no_trees=250, nsim=200000, wishart_df=df, prefix='df,'+str(int(df)), dest_folder=out, sap_sim=sap_sim, sap_ana=sap_ana)
    
def run_posterior_multichain(wishart_df=1000, true_tree_as_identifier=None, result_file='result_mc3.csv', emp_cov_file=None, emp_remove=-1, remove_outgroup=False, make_emp_cov_file=True):
    if true_tree_as_identifier is None:
        true_tree=Rcatalogue_of_trees.tree_good
    else:
        with open(true_tree_as_identifier, 'r') as f:
            s=f.readline().rstrip()
            true_tree=tree_statistics.identifier_to_tree_clean(s)
    if remove_outgroup:
        true_tree=Rtree_operations.remove_outgroup(true_tree)
        true_tree=Rtree_operations.simple_reorder_the_leaves_after_removal_of_s1(true_tree)
    if make_emp_cov_file:
        cov=tree_to_data.get_empirical_matrix(s, factor=0.01, reps=400)
        tree_to_data.emp_cov_to_file(cov, filename=emp_cov_file)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(true_tree)
    no_leaves=Rtree_operations.get_no_leaves(true_tree)
    s_tree=tree_statistics.identifier_to_tree_clean('w.w.a.w.w.a.a.a.w-c.w.c.c.w.w.c.0.w.w.6.3.2-c.w.w.0.w.c.5.w.w-c.w.0.c.3.w.w-c.w.c.2.0-w.c.1-c.0;0.828-0.21-0.197-0.247-0.568-1.06-0.799-1.162-2.632-2.001-0.45-1.048-0.834-0.469-0.191-2.759-0.871-1.896-0.473-0.019-1.236-0.287-0.179-0.981-0.456-0.91-2.114-3.368;0.655-0.506-0.389-0.23')
    print 'no_leaves', no_leaves
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_basic_tree_statistics(tree_statistics.majority_tree, 'majority_tree', output='string'),
               summary.s_variable('add', output='double'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_likelihood(),
               summary.s_prior(),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double_missing') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    if emp_cov_file is not None:
        if emp_remove<0:
            emp_cov=tree_to_data.file_to_emp_cov(emp_cov_file)
        else:
            emp_cov=tree_to_data.file_to_emp_cov(emp_cov_file, emp_remove)
    else:
        emp_cov=None
    print 'emp_cov', emp_cov
    r=simulation_sanity.test_posterior_model_multichain(true_tree, s_tree, [50]*20000, summaries=summaries, thinning_coef=24, wishart_df=wishart_df, result_file=result_file, emp_cov=emp_cov, rescale_empirical_cov=False)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(r)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree)
    
def pack_items(a,b,c):
    res=[]
    for ai in a:
        for bi in b:
            for ci in c:
                res.append([ai,bi,ci])
    return res

def analyse_data_single_chained(filename):
    emp_cov=load_data.read_data(filename, nodes=['French', 'Han', 'Karitiana', 'Sardinian', 'Yoruba'], noss=True)
    print emp_cov
    df=100
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_basic_tree_statistics(tree_statistics.majority_tree, 'majority_tree', output='string'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param'),
               summary.s_variable('rescale_adap_param'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double_missing') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model(None,None, 300000, summaries=summaries, thinning_coef=20, wishart_df= df, emp_cov=emp_cov, no_leaves_true_tree=5)

                
    
def run_posterior_grid(tree_files, alpha, wishart_df):
    
    #true_trees= [tree_generation_laboratory.load_tree(tree_file) for tree_file in tree_files]
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr', output='double_missing'), 
               summary.s_no_admixes(), 
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_likelihood(),
               summary.s_prior()]
    def f(x):
        unsuffixed_filename='.'.join(x.split('.')[:-1])
        true_tree=tree_generation_laboratory.identifier_to_tree_clean_wrapper(tree_generation_laboratory.load_tree(x))
        s_tree=Rtree_operations.create_trivial_tree(Rtree_operations.get_no_leaves(true_tree))
        simulation_sanity.test_posterior_model(true_tree,s_tree, 100, summaries=summaries, thinning_coef=30, 
                                               wishart_df= wishart_df, resimulate_regrafted_branch_length=alpha, 
                                               filename=unsuffixed_filename+'-results.csv')

    from pathos.multiprocessing import Pool
    p=Pool(len(tree_files))
    p.map(f, tree_files)
    
    
    
def run_trivial():
    simulation_sanity.trivial_simulation(8.0,100000)
    print 'sim first'
    summaries=[trivial_mcmc.Trivial_Summary()]
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=2, admixes=0, sim_length=100000, list_of_summaries=summaries)
    analyse_results.full_analysis(summaries,
                  trajectories_for_all_temperatures=False,
                  trajectories_for_all_chains=False,
                  prior_distribution=prior_distribution)
    
if __name__=='__main__':
    #run_posterior_multichain(true_tree_as_identifier='tree2.txt')
    #from argparse import ArgumentParser
    #parser = ArgumentParser(usage='pipeline for admixturebayes', version='0.0.1')
    #parser.add_argument('--df', type=float, default=1000.0, help='degrees of freedom to run under')
    #parser.add_argument('--tree_files', nargs='+', default=['tree.txt'], type=str, help='tree files to be run in parallel')
    #parser.add_argument('--alpha',  default=0, help='the shape parameter of the gamma distribution of the resimulation of the regrafted branch')
    #parser.add_argument('--result_file', type=str, default='', help='file to save results in')
    #options=parser.parse_args()
    
    #run_posterior_grid(['tree2.txt', 'tree3.txt'], 2, 10000)
    import cProfile
     
    print cProfile.run('run_d()')
    #run_analysis_of_proposals()
    #analyse_data_single_chained('example1.treemix_in.gz')
    #run_a()
    #from numpy import array
    #print max_dist(array([[4,3,2]]),array([[9,3,1]]))
    #run_d()
    import sys
    sys.exit()
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser(usage='pipeline for admixturebayes', version='0.0.1')
    parser.add_argument('--df', type=float, default=10.0, help='degrees of freedom to run under')
    parser.add_argument('--sap_simulations', action='store_true',default=False, help='skewed admixture proportion prior in the simulated datasets')
    parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
    parser.add_argument('--true_tree', type=str, default='tree_10_0_3_3.txt', help='file with the true tree to use')
    parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in')
    parser.add_argument('--emp_cov', type= str, default='', help = 'file where the empirical covariance matrix is saved')
    parser.add_argument('--remove_population', type=int, default= -1, help= 'index of the population that should be the root')
    options=parser.parse_args()
    #run_e(options.df, options.result_file, sap_sim=options.sap_simulations, sap_ana=options.sap_analysis)
    
    if options.emp_cov:
        run_posterior_multichain(options.df, options.true_tree, options.result_file, options.emp_cov, options.remove_population, remove_outgroup=True)
    else:
        run_posterior_multichain(options.df, options.true_tree, options.result_file, make_emp_cov_file=False)
