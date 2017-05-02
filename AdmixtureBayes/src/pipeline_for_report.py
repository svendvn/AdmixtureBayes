import simulation_sanity
import analyse_results
import Rcatalogue_of_trees
import Rtree_operations
import summary
import generate_prior_trees
import trivial_mcmc
import tree_statistics


def run_a():
    n=10
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
    #prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, sim_length=nsim, list_of_summaries=list_of_summaries)#, thinning_criteria=max_two)
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
    
def run_d():
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
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param'),
               summary.s_variable('rescale_adap_param'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model(true_tree,s_tree, 200000, summaries=summaries, thinning_coef=3, wishart_df= 1000)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(r)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree)
    
def run_e(df, out):
    simulation_sanity.marginalize_out_data_in_posterior(4, no_trees=250, nsim=200000, wishart_df=df, prefix='df,'+str(int(df)), dest_folder=out)
    
def run_posterior_multichain(wishart_df=1000, true_tree_as_identifier=None, result_file='result_mc3.csv'):
    if true_tree_as_identifier is None:
        true_tree=Rcatalogue_of_trees.tree_good
    else:
        with open(true_tree_as_identifier, 'r') as f:
            s=f.readline().rstrip()
            true_tree=tree_statistics.identifier_to_tree_clean(s)
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
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model_multichain(true_tree,s_tree, [50]*20000, summaries=summaries, thinning_coef=30, wishart_df=wishart_df, admixtures_of_true_tree=None, no_leaves_true_tree=4, result_file=result_file)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(r)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree)
    
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
    #run_posterior_multichain()
    #import sys
    #sys.exit()
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser(usage='pipeline for admixturebayes', version='0.0.1')
    parser.add_argument('--df', type=float, default=100.0, help='degrees of freedom to run under')
    parser.add_argument('--true_tree', type=str, default='', help='file with the true tree to use')
    parser.add_argument('--result_file', type=str, default='', help='file to save results in')
    options=parser.parse_args()
    run_posterior_multichain(options.df, options.true_tree, options.result_file)