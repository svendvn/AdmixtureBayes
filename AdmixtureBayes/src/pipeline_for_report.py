import simulation_sanity
import analyse_results
import Rcatalogue_of_trees
import Rtree_operations
import summary
import generate_prior_trees
import trivial_mcmc


def run_a():
    n=3
    s_tree=Rtree_operations.create_trivial_tree(n)
    summaries=[summary.s_variable('posterior'), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_variable('proposal_type', output='string'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
               
    simulation_sanity.test_prior_model(s_tree, 100000, summaries=summaries, thinning_coef=3)
    def max_two(tree):
        if Rtree_operations.get_number_of_admixes(tree)>2:
            return False
        return True
    list_of_summaries=summaries[2:6]
    nsim=100000
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, sim_length=nsim, list_of_summaries=list_of_summaries)#, thinning_criteria=max_two)
    analyse_results.save_to_csv([tuple(range(nsim))]+[tuple(prior_distribution[summ.name]) for summ in list_of_summaries], list_of_summaries, filename='sim_prior.csv', origin_layer=None)
    analyse_results.generate_summary_csv(summaries)
    #analyse_results.full_analysis(summaries,
    #              trajectories_for_all_temperatures=False,
    #              trajectories_for_all_chains=False,
    #              prior_distribution=prior_distribution)
    
def run_b():
    n=4
    s_tree=Rtree_operations.create_trivial_tree(n)
    summaries=[summary.s_variable('posterior'), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_variable('proposal_type', output='string')]+[summary.s_variable(s,output='double_missing') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    simulation_sanity.test_prior_model_no_admixes(s_tree, 100000, summaries=summaries, thinning_coef=3)
    list_of_summaries=summaries[2:6]
    nsim=20000
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, admixes=0, sim_length=nsim, list_of_summaries=list_of_summaries)
    analyse_results.save_to_csv([tuple(range(nsim))]+[tuple(prior_distribution[summ.name]) for summ in list_of_summaries], list_of_summaries, filename='sim_prior.csv', origin_layer=None)
    analyse_results.generate_summary_csv(summaries)
    #analyse_results.full_analysis(summaries,
    #              trajectories_for_all_temperatures=False,
    #              trajectories_for_all_chains=False,
    #              prior_distribution=prior_distribution)
    
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
    run_a()