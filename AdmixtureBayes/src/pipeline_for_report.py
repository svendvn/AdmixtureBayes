import simulation_sanity
import analyse_results
import Rcatalogue_of_trees
import Rtree_operations
import summary
import generate_prior_trees
import trivial_mcmc
import tree_statistics
import tree_generation_laboratory
import load_data


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
               
    #simulation_sanity.test_prior_model(s_tree, 50000, summaries=summaries, thinning_coef=3)
    def max_two(tree):
        if Rtree_operations.get_number_of_admixes(tree)>2:
            return False
        return True
    list_of_summaries=summaries[2:10]
    nsim=100000
    prior_distribution=generate_prior_trees.get_distribution_under_prior(leaves=n, sim_length=nsim, list_of_summaries=list_of_summaries, skewed_admixture_prior=True)#, thinning_criteria=max_two)
    analyse_results.save_to_csv([tuple(range(nsim))]+[tuple(prior_distribution[summ.name]) for summ in list_of_summaries], list_of_summaries, filename='sim_prior2.csv', origin_layer=None)
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
    true_tree=generate_prior_trees.generate_phylogeny(8,2)
    #tree_statistics.identifier_to_tree_clean('w.c.w.1-w.a.w-c.0.w.w-c.w.0-c.0;0.578-1.783-0.221-1.462-0.44-0.295-0.142-0.255-0.352;0.373')
    s_tree=Rtree_operations.create_trivial_tree(8)
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
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model(true_tree,s_tree, 3000, summaries=summaries, thinning_coef=20, wishart_df= 1000, resimulate_regrafted_branch_length=False, admixtures_of_true_tree=2, no_leaves_true_tree=8)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(r)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree)
    
def run_e(df, out, sap_sim, sap_ana):
    simulation_sanity.marginalize_out_data_in_posterior(4, no_trees=250, nsim=200000, wishart_df=df, prefix='df,'+str(int(df)), dest_folder=out, sap_sim=sap_sim, sap_ana=sap_ana)
    
def run_posterior_multichain(wishart_df=1000, true_tree_as_identifier=None, result_file='result_mc3.csv'):
    if true_tree_as_identifier is None:
        true_tree=Rcatalogue_of_trees.tree_good
    else:
        with open(true_tree_as_identifier, 'r') as f:
            s=f.readline().rstrip()
            true_tree=tree_statistics.identifier_to_tree_clean(s)
    print 'true_tree', tree_statistics.unique_identifier_and_branch_lengths(true_tree)
    no_leaves=Rtree_operations.get_no_leaves(true_tree)
    s_tree=Rtree_operations.create_trivial_tree(no_leaves)
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
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_likelihood(),
               summary.s_prior(),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model_multichain(true_tree, s_tree, [50]*2000, summaries=summaries, thinning_coef=24, wishart_df=wishart_df, result_file=result_file)
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
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s,output='double') for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    r=simulation_sanity.test_posterior_model(None,None, 300000, summaries=summaries, thinning_coef=20, wishart_df= df, emp_cov=emp_cov, no_leaves_true_tree=5)

                
    
def run_posterior_grid(tree_files, alpha, wishart_df):
    
    #true_trees= [tree_generation_laboratory.load_tree(tree_file) for tree_file in tree_files]
    summaries=[summary.s_posterior(), 
               summary.s_variable('mhr'), 
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
    #import cProfile
     
    #print cProfile.run('run_d()')
    #run_d()
    #analyse_data_single_chained('example1.treemix_in.gz')
    #import sys
    #sys.exit()
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser(usage='pipeline for admixturebayes', version='0.0.1')
    parser.add_argument('--df', type=float, default=1000.0, help='degrees of freedom to run under')
    parser.add_argument('--sap_simulations', action='store_true',default=False, help='skewed admixture proportion prior in the simulated datasets')
    parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
    parser.add_argument('--true_tree', type=str, default='tree3.txt', help='file with the true tree to use')
    parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in')
    options=parser.parse_args()
    #run_e(options.df, options.result_file, sap_sim=options.sap_simulations, sap_ana=options.sap_analysis)
    
    run_posterior_multichain(options.df, options.true_tree, options.result_file)