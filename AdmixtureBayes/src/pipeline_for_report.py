import simulation_sanity
import analyse_results
import Rcatalogue_of_trees
import Rtree_operations
import summary
import generate_prior_trees

def run_a():
    s_tree=Rtree_operations.create_trivial_tree(4)
    simulation_sanity.test_prior_model(s_tree, 100000)
    summaries=[summary.s_variable('posterior'), summary.s_variable('mhr'), summary.s_no_admixes()]
    distr=generate_prior_trees.get_distribution_under_prior(leaves=4, sim_length=1000, list_of_summaries=[summaries[2]])
    prior_distribution={summaries[2].name: distr[0]}
    analyse_results.full_analysis(summaries,
                  trajectories_for_all_temperatures=False,
                  trajectories_for_all_chains=False,
                  prior_distribution=prior_distribution)
    
if __name__=='__main__':
    run_a()