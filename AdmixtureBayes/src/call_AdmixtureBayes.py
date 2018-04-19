from argparse import ArgumentParser

from construct_proposal_choices import make_proposal
from construct_starting_trees_choices import get_starting_trees
from construct_covariance_choices import get_covariance
from construct_nodes_choices import get_nodes
from construct_summary_choices import get_summary_scheme
from construct_filter_choices import make_filter
from temperature_scheme import fixed_geometrical
from analyse_results import save_permuts_to_csv, get_permut_filename
from posterior import posterior_class
from MCMCMC import MCMCMC
from wishart_distribution_estimation import estimate_degrees_of_freedom
from MCMC import basic_chain
from stop_criteria import stop_criteria
from one_evaluation import one_evaluation
from tree_to_data import emp_cov_to_file, file_to_emp_cov


import os 
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

#famous tree:
#w.w.w.w.w.w.a.a.w-c.w.c.c.w.c.5.0.w.3.2-c.w.w.0.c.4.w-c.w.0.c.3-w.c.1-c.0;0.07-0.974-1.016-0.089-0.81-0.086-1.499-0.052-1.199-2.86-0.403-0.468-0.469-1.348-1.302-1.832-0.288-0.18-0.45-0.922-2.925-3.403;0.388-0.485

parser = ArgumentParser(usage='pipeline for Admixturebayes', version='1.0.0')

#input/output options
parser.add_argument('--input_file', type=str, default='(4,2)', help='the input file of the pipeline. Its type should match the first argument of covariance_pipeline. 6= treemix file, 7-9=covariance file')
parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in. The prefix will not be prepended the result_file.')
parser.add_argument('--prefix', type=str, default='sletmig/', help= 'this directory will be the beginning of every temporary file created in the covariance pipeline and in the estimation of the degrees of freedom in the wishart distribution.')

#special run options
parser.add_argument('--profile', action='store_true', default=False, help="this will embed the MCMC part in a profiler")
parser.add_argument('--treemix_instead', action= 'store_true', default=False, help='this will call treemix instead of AdmixtureBayes')
parser.add_argument('--treemix_also', action='store_true', default=False, help='this will call treemix in addition to AdmixtureBayes')
parser.add_argument('--likelihood_treemix', action='store_true', default=False, help='this will use the likelihood from treemix instead of the wishart distribution.')
parser.add_argument('--evaluate_likelihood', action='store_true', default=False, help='this will evaluate the likelihood in the starting tree and then stop, writing just a single file with three values, prior, likelihood and posterior.')
parser.add_argument('--evaluate_bootstrap_likelihoods', action='store_true', default=False, help='If evaluate likelihood is turned on this will calculate the likelihood of all bootstrapped covariances(if bootstrapping is also turned on)')
parser.add_argument('--stop_evaluations', action='store_true', default=False, help='This will stop the analysis after the data preparation')

#treemix arguments
parser.add_argument('--treemix_reps', type=int, default=1, help='the number of repititions of the treemix call. Only used when treemix_instead or treemix_also')
parser.add_argument('--treemix_no_admixtures', type=int, nargs='+', default=[0,1,2,3], help='the number of admixture events in treemixrun. Only used when treemix_instead or treemix_also')
parser.add_argument('--treemix_processes', type=int, default=1, help='the number of parallel processes to run treemix over.')
parser.add_argument('--alternative_treemix_infile', type=str, default='', help='By default the program will use the treemix file generated in the covariance pipeline (or go looking for the file that would have been made if 6 was part of the pipeline). This will override that')
#parser.add_argument('--treemix_file', type=str, default='', help= 'the filename of the intermediate step that contains the ms output.')
parser.add_argument('--treemix_output_prefix', type=str, default='', help= 'the filename prefix of all the treemix output files. Each file will get the suffix k.txt where k is the number of admixture events.')

#covariance matrix options
parser.add_argument('--covariance_pipeline', nargs='+', type=int, default=[2,3,4,5,6,7,8,9], help='skewed admixture proportion prior in the simulated datasets')
parser.add_argument('--outgroup_name', type=str, default='out', help='the name of the outgroup that should be added to a simulated dataset.')
parser.add_argument('--reduce_node', type=str, default='out', help='the name of the population that should be made root.')
parser.add_argument('--nodes', type=str, nargs='+', default=[''], help= 'list of nodes of the populations or the filename of a file where the first line contains all population names. If unspecified the first line of the input_file will be used. If no input file is found, there will be used standard s1,..,sn.')
parser.add_argument('--wishart_noise', action='store_true', default=False)

#empirical_covariance matrix. Base estimation
parser.add_argument('--arcsin', action='store_true', default=False)
parser.add_argument('--cov_estimation', choices=['None', 'Jade','outgroup_sum', 'outgroup_product', 'average_sum', 'average_product','Jade-o', 'EM'], default='average_sum', help='this is the way of estimating the empirical covariance matrix.')
parser.add_argument('--bias_c_weight', choices=['default','None','outgroup_sum', 'outgroup_product', 'average_sum', 'average_product'], default='default', help='from cov_weight with bias correction unweighted there are some obvious choices for weighing the bias correction, so here they are: None=None, Jade=average_sum, Jade-o=outgroup_sum, average_sum=average_sum, average_product=average_product, outgroup_sum=outgroup_sum, outgroup_product=outgroup_product')
parser.add_argument('--Jade_cutoff', type=float, default=1e-5, help='this will remove SNPs of low diversity in either the Jade or the Jade-o scheme.')
parser.add_argument('--variance_correction', default='unbiased', choices=['None', 'unbiased','mle'], help= 'The type of adjustment used on the empirical covariance.')
parser.add_argument('--variance_correction_input_file', default='', type=str, help='if the variance correction is saved in a file (with numpy.savetxt format of a 2 dimensional numpy array) it can be loaded in with this command')
parser.add_argument('--add_variance_correction_to_graph', default=True, action='store_true', help='If on, the variance correction will be added to the covariance matrix of the graph and not subtracted from the empirical covariance matrix.')
parser.add_argument('--indirect_correction', default=False, action='store_true', help='the bias in the covariance is (possibly again) corrected for by indirect estimation.')
parser.add_argument('--indirect_its', type=int, default=100, help='For how many iterations should the indirect optimization procedure be run. Only applicable if indirect_correction is True')
parser.add_argument('--indirect_simulation_factor', type=int, default=1, help='How much more data than provided should be simulated in the indirect correction procedure. Only applicable if indirect_correction is True')
parser.add_argument('--EM_maxits', type=int, default=100, help='The maximum number of iterations of the EM algorithm. There is another stopping criteria that may stop it before.')
parser.add_argument('--EM_alpha', type=float, default=1.0, help='The EM algorithm assumes that the allele frequencies in the outgroup are known. In fact it is estimated with: if alpha=1.0: the empirical allele frequencies of the outgroup, alpha=0.0: the average empirical allele frequency in all the other populations.It can also be chosen as something in between. This estimator is biased because the outgroup allele frequencies are not known and because of extra normal distribution assumptions')
parser.add_argument('--no_repeats_of_cov_est', type=int, default=1, help='The number of times the simulation procedure should be run.')
parser.add_argument('--indirect_randomize_seed', action='store_true', default=False, help='This will make indirect estimation (if indirect_correction) use different seeds, slowing down and making maximization more troublesome yet being more correct.')
parser.add_argument('--initial_Sigma', choices=['default','random', 'start'], default='default', help='This means that ')
#empirical_covariance matrix. Post estimation
parser.add_argument('--scale_goal', choices=['min','max'], default='max', help='If 9 is included in the pipeline, this is what there will be scaled to.')

#prior options
parser.add_argument('--p', type=float, default=0.5, help= 'the geometrical parameter in the prior. The formula is p**x(1-p)')
parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
parser.add_argument('--uniform_prior', action='store_true', default=False, help='If applied a uniform prior will be used on the different topologies.')
parser.add_argument('--no_add', action='store_true', default=False, help='this will remove the add contribution')

#filter arguments, if 6 and higher is in the covariance pipeline, there will be applied a filter to the SNPs
parser.add_argument('--filter_type', choices=['snp','none', 'outgroup_other','outgroup','all_pops'], default='snp', help='This will apply a filter to positions based on their value.')
parser.add_argument('--filter_on_simulated', choices=['same', 'none', 'outgroup_other', 'outgroup', 'snp', 'all_pops'], default='same', help='In indirect inference, whole datasets are simulated under ')

#degrees of freedom arguments
parser.add_argument('--estimate_bootstrap_df', action='store_true', default=False, help= 'if declared, the program will estimate the degrees of freedom in the wishart distribution with a bootstrap sample.')
parser.add_argument('--df_file', type=str, default='', help='file with the number of degrees of freedom in the wishart distribution OR the file of entrywise variances of the covariance matrix.')
parser.add_argument('--wishart_df', type=float, default=1000.0, help='degrees of freedom to run under if bootstrap-mle of this number is declined (and treemix_likelihood is not specified).')
parser.add_argument('--save_df_file', type=str, default='DF.txt', help='the prefix is put before this string and the degrees of freedom is saved to this file.')
parser.add_argument('--bootstrap_blocksize', type=int, default=1000, help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')
parser.add_argument('--no_bootstrap_samples', type=int, default=100, help='the number of bootstrap samples to make to estimate the degrees of freedom in the wishart distribution.')
parser.add_argument('--df_treemix_adjust_to_wishart', action='store_true', default=False, help='This will, if likelihood_treemix is flagged and df_file is a wishart-df, choose a variance matrix that gives a normal distribution with the same mode-likelihood-value as if no likelihood_treemix had been switched on.')
parser.add_argument('--save_bootstrap_covariances', type=str, default='', help='if provided the bootstrapped covariance matrices will be saved to numbered files starting with {prefix}+_+{save_covariances}+{num}+.txt')
parser.add_argument('--bootstrap_type_of_estimation', choices=['mle_opt','var_opt'], default='var_opt', help='This is the way the bootstrap wishart estimate is estimated.')
parser.add_argument('--load_bootstrapped_covariances', type=str, default=[], nargs='+', help='if supplied, this will load covariance matrices from the specified files instead of simulating new ones.')

#proposal frequency options
parser.add_argument('--deladmix', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--addadmix', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--regraft', type=float, default=0, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_add', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_admix', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_admix_correction', type=float, default=0, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_constrained', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_marginally', type=float, default=0, help='this states the frequency of the proposal type')
parser.add_argument('--sliding_regraft', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--sliding_rescale', type=float, default=0, help='this states the frequency of the proposal type')

#other proposal options
parser.add_argument('--cancel_preserve_root_distance', default=False, action='store_true', help="if applied there will not be made correction for root distance when adding and deleting admixtures")

#mc3 arguments
parser.add_argument('--starting_trees', type=str, nargs='+', default=[], help='filenames of trees to start in. If empty, the trees will either be simulated with the flag --random_start or the so-called trivial tree')
parser.add_argument('--starting_adds', type=str, nargs='+', default=[], help="filename of the adds to use on the starting trees.")
parser.add_argument('--start', choices=['trivial','random', 'perfect'], default='trivial', help='Where to start the chain - works only if starting trees are not specified.')
parser.add_argument('--starting_tree_scaling', choices=['None','empirical_trace', 'starting_tree_trace','scalar'], default='None', type=str, help='The starting tree can be scaled as the covariance (as_covariance) or as the p')
parser.add_argument('--starting_tree_use_scale_tree_factor', default=False, action='store_true', help='this will scale the tree with the specified scale_tree_factor.')

#tree simulation
parser.add_argument('--p_sim', type=float, default=.5, help='the parameter of the geometric distribution in the distribution to simulate the true tree from.')
parser.add_argument('--popsize', type=int, default=20, help='the number of genomes sampled from each population.')
parser.add_argument('--nreps', type=int, default=50, help='How many pieces of size 500 kb should be simualted')
parser.add_argument('--scale_tree_factor', type=float, default=0.02, help='The scaling factor of the simulated trees to make them less vulnerable to the fixation effect.')
parser.add_argument('--skewed_admixture_prior_sim', default=False, action='store_true', help='the prior tree is simulated with an uneven prior on the admixture proportions')
parser.add_argument('--time_adjusted_tree', default=False, action='store_true', help='this will modify the simulated tree such that all drift lengths from root to leaf are the same')
parser.add_argument('--sadmix_tree', default=False, action='store_true', help='this will simulate trees where all admixture events are important in the sense that they expand the space of possible covariance matrices.')

#covariance simulation
parser.add_argument('--favorable_init_brownian', default=False, action='store_true', help='This will start the brownian motion(only if 21 in workflow) between 0.4 and 0.6')
parser.add_argument('--unbounded_brownian', default=False, action='store_true', help='This will start the brownian motion(only if 21 in workflow) between 0.4 and 0.6')
parser.add_argument('--filter_on_outgroup', default=False, action='store_true', help='If applied (and 23 in the pipeline) SNPs that are not polymorphic in the outgroup are removed. If not, the default is that polymorphic in no population are removed. ')

#chain data collection
parser.add_argument('--summary_majority_tree', action='store_true', default=False, help='this will calculate the majority (newick) tree based on the sampled tree')
parser.add_argument('--summary_acceptance_rate', action='store_true', default=True, help='This will calculate and store summaries related to the acceptance rate')
parser.add_argument('--summary_admixture_proportion_string', action='store_true', default=True, help='this will save a string in each step indicating names and values of all admixture proportions')

#MCMCMC setup
parser.add_argument('--MCMC_chains', type=int, default=8, help='The number of chains to run the MCMCMC with.')
parser.add_argument('--n', type=int, default=200, help='the number of MCMCMC flips throughout the chain.')
parser.add_argument('--m', type=int, default=50, help='the number of MCMC steps before the chain is ')
parser.add_argument('--max_temp', type=float, default=100, help='the maximum temperature used in the MCMCMC.')
parser.add_argument('--thinning_coef', type=int, default=40, help='the number of iterations between each data recording point.')
parser.add_argument('--store_permuts', action='store_true', default=False, help='If applied, the permutations from the MCMCMC flips are recorded in a file with a similar filename to the result_file')
parser.add_argument('--stop_criteria', action='store_true', default=True, help='If applied the MCMCMC will stop when the coldest chain has an effective sample size at ')
parser.add_argument('--stop_criteria_frequency', type=int, default=200000, help='This tells the frequency of checking for when the stop criteria are checked (if the stop_criteria flag is turned on)')
parser.add_argument('--stop_criteria_topological', default=False, action='store_true', help='this will add a topological stop criteria that should also be fulfilled.')

options=parser.parse_args()





mp=make_proposal(options)

before_added_outgroup, full_nodes, reduced_nodes=get_nodes(options.nodes, options.input_file, options.outgroup_name, options.reduce_node)
print 'before_nodes', before_added_outgroup
print 'full_nodes', full_nodes
print 'reduced_nodes', reduced_nodes




if options.prefix and options.prefix[-1]!='_':
    prefix=options.prefix+'_'
else:
    prefix=options.prefix

if options.alternative_treemix_infile:
    treemix_file=options.alternative_treemix_infile
    treemix_in_file=options.alternative_treemix_infile
else:
    treemix_file=prefix+"treemix_in.txt"
    treemix_in_file=treemix_file+'.gz'
    
#treemix_file=prefix+"treemix_in.txt" 
if options.treemix_output_prefix:
    treemix_out_files=prefix+options.treemix_output_prefix
else:
    treemix_out_files=prefix+'treemix_out'

if options.starting_trees:
    preliminary_starting_trees=get_starting_trees(options.starting_trees, 
                                      options.MCMC_chains, 
                                      adds=options.starting_adds,
                                      nodes=reduced_nodes, 
                                      pipeline=options.covariance_pipeline,
                                      multiplier=None,
                                      scale_tree_factor=options.scale_tree_factor,
                                      start=options.start, 
                                      prefix=prefix,
                                      starting_tree_scaling=options.starting_tree_scaling,
                                      starting_tree_use_scale_tree_factor=options.starting_tree_use_scale_tree_factor,
                                      scale_goal=options.scale_goal)
else:
    preliminary_starting_trees=[None]
    assert options.initial_Sigma!='start', 'to make the filter start somewhere specific it should also be specified specifically'

locus_filter=make_filter(options.filter_type,
                         outgroup_name=options.reduce_node,
                         covariance_pipeline=options.covariance_pipeline)
if options.filter_on_simulated=='same':
    locus_filter_on_simulated=make_filter(options.filter_type)
else:
    locus_filter_on_simulated=make_filter(options.filter_on_simulated)
    
    


estimator_arguments=dict(reducer=options.reduce_node,
                         variance_correction=options.variance_correction,
                         method_of_weighing_alleles=options.cov_estimation,
                         arcsin_transform=options.arcsin,
                         jade_cutoff=options.Jade_cutoff,
                         bias_c_weight=options.bias_c_weight,
                         nodes=full_nodes,
                         indirect_correction=options.indirect_correction,
                         Indirect_its=options.indirect_its,
                         Indirect_multiplier_s=options.indirect_simulation_factor,
                         EM_maxits=options.EM_maxits,
                         EM_alpha=options.EM_alpha,
                         no_repeats_of_cov_est=options.no_repeats_of_cov_est,
                         Simulator_fixed_seed=not options.indirect_randomize_seed,
                         initial_Sigma_generator={options.initial_Sigma:(preliminary_starting_trees[0], reduced_nodes)},
                         locus_filter_on_simulated=locus_filter_on_simulated,
                         add_variance_correction_to_graph=options.add_variance_correction_to_graph,
                         save_variance_correction=True,
                         prefix=prefix)

covariance=get_covariance(options.covariance_pipeline, 
                          options.input_file, 
                          full_nodes=full_nodes, 
                          skewed_admixture_prior_sim=options.skewed_admixture_prior_sim, 
                          p=options.p_sim, 
                          outgroup_name=options.outgroup_name,
                          reduce_covariance_node=options.reduce_node,
                          sample_per_pop=options.popsize,
                          nreps=options.nreps,
                          treemix_file=treemix_file,
                          scale_tree_factor=options.scale_tree_factor,
                          prefix=prefix,
                          t_adjust_tree=options.time_adjusted_tree,
                          sadmix=options.sadmix_tree,
                          add_wishart_noise_to_covariance=options.wishart_noise,
                          df_of_wishart_noise_to_covariance=options.wishart_df,
                          scale_goal=options.scale_goal,
                          favorable_init_brownian=options.favorable_init_brownian,
                          unbounded_brownian=options.unbounded_brownian,
                          filter_on_outgroup=options.filter_on_outgroup,
                          locus_filter=locus_filter,
                          estimator_arguments=estimator_arguments
                          )

if options.treemix_instead or options.treemix_also:

    dir = os.path.dirname(__file__)
    program=os.path.join(dir,'treemixrunner.py')
    calls_to_treemix=[['python',program, 
                       '-p', str(options.treemix_processes), 
                       '-n', str(options.treemix_reps), 
                       '-i', treemix_in_file,
                       '-o', treemix_out_files+str(k),
                       '-m', str(k)] for k in options.treemix_no_admixtures]
    from subprocess import call
    for c in calls_to_treemix:
        call(c)
    
    if not options.treemix_also:
        import sys
        sys.exit()



if options.estimate_bootstrap_df:
    #assert 6 in options.covariance_pipeline, 'Can not estimate the degrees of freedom without SNP data.'
    #reduce_also= (8 in options.covariance_pipeline)
    estimator_arguments['save_variance_correction']=False
    if options.likelihood_treemix:
        summarization='var'
    else:
        summarization=options.bootstrap_type_of_estimation
    df, boot_covs=estimate_degrees_of_freedom(treemix_in_file, 
                                           bootstrap_blocksize=options.bootstrap_blocksize, 
                                           no_bootstrap_samples=options.no_bootstrap_samples,
                                           summarization=summarization,
                                           cores=options.MCMC_chains,
                                           save_covs=options.save_bootstrap_covariances,
                                           prefix=prefix,
                                           est=estimator_arguments, locus_filter=locus_filter,
                                           load_bootstrapped_covariances=options.load_bootstrapped_covariances)
elif options.df_file:
    if options.likelihood_treemix and not options.df_treemix_adjust_to_wishart:
        df=file_to_emp_cov(options.df_file, nodes=reduced_nodes)
    else:
        with open(options.df_file, 'r') as f:
            df=float(f.readline().rstrip())
    boot_covs=[]
else:
    df=options.wishart_df
    boot_covs=[]




if 9 not in options.covariance_pipeline:
    multiplier=None
    covariance=(covariance, multiplier)
else:
    multiplier=covariance[1]


starting_trees=get_starting_trees(options.starting_trees, 
                                  options.MCMC_chains, 
                                  adds=options.starting_adds,
                                  nodes=reduced_nodes, 
                                  pipeline=options.covariance_pipeline,
                                  multiplier=multiplier,
                                  scale_tree_factor=options.scale_tree_factor,
                                  start=options.start, 
                                  prefix=prefix,
                                  starting_tree_scaling=options.starting_tree_scaling,
                                  starting_tree_use_scale_tree_factor=options.starting_tree_use_scale_tree_factor,
                                  scale_goal=options.scale_goal)

# if not options.starting_trees:
#     no_pops=len(reduced_nodes)
#     starting_trees=map(str, [no_pops]*options.MCMC_chains)
# else:
#     starting_trees=options.starting_trees

    


if (not options.likelihood_treemix) or options.df_treemix_adjust_to_wishart:
    with open(prefix+options.save_df_file, 'w') as f:
        f.write(str(df))
else:
    emp_cov_to_file(df, prefix+options.save_df_file, nodes=reduced_nodes)
    

summary_verbose_scheme, summaries=get_summary_scheme(majority_tree=options.summary_majority_tree, 
                                          light_newick_tree_summaries=options.stop_criteria_topological,
                                          full_tree=True, #can not think of a moment where you don't want this.
                                          proposals=mp[0], 
                                          acceptance_rate_information=options.summary_acceptance_rate,
                                          admixture_proportion_string=options.summary_admixture_proportion_string,
                                          no_chains=options.MCMC_chains)

sim_lengths=[options.m]*options.n


    
# print 'starting_trees', starting_trees
# print 'posterior', posterior
# print 'summaries',summaries
# print 'temperature_scheme', fixed_geometrical(options.max_temp,options.MCMC_chains)
# print 'summary_verbose_scheme',summary_verbose_scheme
# #print 'sim_lengths',sim_lengths
# print 'int(options.thinning_coef)',int(options.thinning_coef)
# print 'mp',mp
# print 'options.MCMC_chains',options.MCMC_chains
# print 'multiplier', multiplier
# print 'result_file', options.result_file
# print 'options.store_permuts', options.store_permuts

if options.stop_criteria:
    sc=stop_criteria(frequency=options.stop_criteria_frequency, outfile=prefix+'stop_criteria.txt', topological=options.stop_criteria_topological)
else:
    sc=None
    
if options.stop_evaluations:
    import sys
    sys.exit()
    
posterior= posterior_class(emp_cov=covariance[0], 
                       M=df, 
                       p=options.p, 
                       use_skewed_distr=options.sap_analysis, 
                       multiplier=covariance[1], 
                       nodes=reduced_nodes, 
                       use_uniform_prior=options.uniform_prior, 
                       treemix=options.likelihood_treemix,
                       add_variance_correction_to_graph=options.add_variance_correction_to_graph,
                       prefix=prefix,
                       variance_correction_file=options.variance_correction_input_file)
    

def multi_chain_run():
    res=MCMCMC(starting_trees=starting_trees, 
           posterior_function= posterior,
           summaries=summaries, 
           temperature_scheme=fixed_geometrical(options.max_temp,options.MCMC_chains), 
           printing_schemes=summary_verbose_scheme, 
           iteration_scheme=sim_lengths, 
           overall_thinnings=int(options.thinning_coef), 
           proposal_scheme= mp, 
           cores=options.MCMC_chains, 
           no_chains=options.MCMC_chains,
           multiplier=multiplier,
           result_file=options.result_file,
           store_permuts=options.store_permuts, 
           stop_criteria=sc)
    
def single_chain_run():
    basic_chain(start_x= starting_trees[0],
                summaries=summaries,
                posterior_function=posterior,
                proposal=mp[0],
                post=None,
                N=sum(sim_lengths),
                sample_verbose_scheme=summary_verbose_scheme[0],
                overall_thinning=int(options.thinning_coef),
                i_start_from=0,
                temperature=1.0,
                proposal_update=None,
                multiplier=multiplier,
                appending_result_file=options.result_file,
                appending_result_frequency=sim_lengths[0])
    
def single_evaluation_run(alts=[]):
    one_evaluation(starting_trees[0], 
                   posterior, 
                   options.result_file,
                   alts)
    
if options.evaluate_likelihood:
    if options.evaluate_bootstrap_likelihoods and boot_covs:
        single_evaluation_run(boot_covs)
    else:
        single_evaluation_run()
    with open(options.result_file, 'r') as f:
        for r in f.readlines():
            print r
else:
    if options.profile:
        import cProfile
        cProfile.run('single_chain_run()')
    elif options.MCMC_chains==1:
        single_chain_run()
    else:
        res=multi_chain_run()
        if options.store_permuts:
            save_permuts_to_csv(res, get_permut_filename(options.result_file))


        
