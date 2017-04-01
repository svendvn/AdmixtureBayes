import simulation_sanity
import analyse_results
import Rcatalogue_of_trees
import Rtree_operations
import summary
import generate_prior_trees

def run_a():
    n=3
    s_tree=Rtree_operations.create_trivial_tree(n)
    summaries=[summary.s_variable('posterior'), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier()]
    simulation_sanity.test_prior_model(s_tree, 100000, summaries=summaries)
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, sim_length=1000, list_of_summaries=summaries[2:])
    analyse_results.full_analysis(summaries,
                  trajectories_for_all_temperatures=False,
                  trajectories_for_all_chains=False,
                  prior_distribution=prior_distribution)
    
def run_b():
    n=4
    s_tree=Rtree_operations.create_trivial_tree(n)
    summaries=[summary.s_variable('posterior'), 
               summary.s_variable('mhr'), 
               summary.s_tree_identifier()]
    simulation_sanity.test_prior_model_no_admixes(s_tree, 100000, summaries=summaries)
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, admixes=0, sim_length=1000, list_of_summaries=summaries[2:])
    analyse_results.full_analysis(summaries,
                  trajectories_for_all_temperatures=False,
                  trajectories_for_all_chains=False,
                  prior_distribution=prior_distribution)
    
if __name__=='__main__':
    run_a()