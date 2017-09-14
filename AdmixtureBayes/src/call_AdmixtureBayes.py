from argparse import ArgumentParser
from meta_proposal import simple_adaptive_proposal
from construct_starting_trees_choices import get_starting_trees
from construct_covariance_choices import get_covariance
from construct_nodes_choices import get_nodes
from construct_summary_choices import get_summary_scheme
from temperature_scheme import fixed_geometrical
from analyse_results import save_permuts_to_csv, get_permut_filename
from posterior import initialize_posterior
from MCMCMC import MCMCMC
from wishart_distribution_estimation import estimate_degrees_of_freedom
from MCMC import basic_chain
from stop_criteria import stop_criteria

import os 
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

#famous tree:
#w.w.w.w.w.w.a.a.w-c.w.c.c.w.c.5.0.w.3.2-c.w.w.0.c.4.w-c.w.0.c.3-w.c.1-c.0;0.07-0.974-1.016-0.089-0.81-0.086-1.499-0.052-1.199-2.86-0.403-0.468-0.469-1.348-1.302-1.832-0.288-0.18-0.45-0.922-2.925-3.403;0.388-0.485

parser = ArgumentParser(usage='pipeline for Admixturebayes', version='1.0.0')

#overall options
parser.add_argument('--input_file', type=str, default='(6,0)', help='the input file of the pipeline. Its type should match the first argument of covariance_pipeline. 6= treemix file, 7-9=covariance file')
parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in. The prefix will not be prepended the result_file.')
parser.add_argument('--prefix', type=str, default='sletmig/', help= 'this directory will be the beginning of every temporary file created in the covariance pipeline and in the estimation of the degrees of freedom in the wishart distribution.')
parser.add_argument('--profile', action='store_true', default=False, help="this will embed the MCMC part in a profiler")
parser.add_argument('--treemix_instead', action= 'store_true', default=False, help='this will call treemix instead of AdmixtureBayes')
parser.add_argument('--treemix_also', action='store_true', default=False, help='this will call treemix in addition to AdmixtureBayes')

#treemix arguments
parser.add_argument('--treemix_reps', type=int, default=1, help='the number of repititions of the treemix call. Only used when treemix_instead or treemix_also')
parser.add_argument('--treemix_no_admixtures', type=int, nargs='+', default=[0,1,2,3], help='the number of admixture events in treemixrun. Only used when treemix_instead or treemix_also')
parser.add_argument('--treemix_processes', type=int, default=1, help='the number of parallel processes to run treemix over.')
parser.add_argument('--alternative_treemix_infile', type=str, default='', help='By default the program will use the treemix file generated in the covariance pipeline (or go looking for the file that would have been made if 6 was part of the pipeline). This will override that')


#covariance matrix options
parser.add_argument('--covariance_pipeline', nargs='+', type=int, default=[2,3,4,5,6,7,8,9], help='skewed admixture proportion prior in the simulated datasets')
parser.add_argument('--outgroup_name', type=str, default='out', help='the name of the outgroup that should be added to a simulated dataset.')
parser.add_argument('--reduce_node', type=str, default='out', help='the name of the population that should be made root.')
parser.add_argument('--nodes', type=str, nargs='+', default=[''], help= 'list of nodes of the populations or the filename of a file where the first line contains all population names. If unspecified the first line of the input_file will be used. If no input file is found, there will be used standard s1,..,sn.')
parser.add_argument('--wishart_noise', action='store_true', default=False)

#prior options
parser.add_argument('--p', type=float, default=0.5, help= 'the geometrical parameter in the prior. The formula is p**x(1-p)')
parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
parser.add_argument('--uniform_prior', action='store_true', default=False, help='If applied a uniform prior will be used on the different topologies.')
parser.add_argument('--no_add', action='store_true', default=False, help='this will remove the add contribution')

#degrees of freedom arguments
parser.add_argument('--estimate_bootstrap_df', action='store_true', default=False, help= 'if declared, the program will estimate the degrees of freedom in the wishart distribution with a bootstrap sample.')
parser.add_argument('--wishart_df', type=float, default=1000.0, help='degrees of freedom to run under if bootstrap-mle of this number is declined.')
parser.add_argument('--bootstrap_blocksize', type=int, default=1000, help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')
parser.add_argument('--no_bootstrap_samples', type=int, default=100, help='the number of bootstrap samples to make to estimate the degrees of freedom in the wishart distribution.')

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
parser.add_argument('--MCMC_chains', type=int, default=7, help='The number of chains to run the MCMCMC with.')
parser.add_argument('--starting_trees', type=str, nargs='+', default=[], help='filenames of trees to start in. If empty the trees will either be simulated with the flag --random_start or the so-called trivial tree')
parser.add_argument('--starting_add', type=str, default='', help="filename of file with the optimal add")
parser.add_argument('--random_start', action='store_true', default=False, help='If supplied, the starting trees will be simulated from the prior (unless the starting trees are specified)')

#tree simulation
parser.add_argument('--p_sim', type=float, default=.5, help='the parameter of the geometric distribution in the distribution to simulate the true tree from.')
parser.add_argument('--popsize', type=int, default=20, help='the number of genomes sampled from each population.')
parser.add_argument('--nreps', type=int, default=500, help='How many pieces of size 500 kb should be simualted')
parser.add_argument('--treemix_file', type=str, default='', help= 'the filename of the intermediate step that contains the ms output.')
parser.add_argument('--ms_variance_correction', default=False, action='store_true', help= 'Should the empirical covariance matrix be adjusted for finite sample size.')
parser.add_argument('--scale_tree_factor', type=float, default=0.02, help='The scaling factor of the simulated trees to make them less vulnerable to the fixation effect.')
parser.add_argument('--skewed_admixture_prior_sim', default=False, action='store_true', help='the prior tree is simulated with an uneven prior on the admixture proportions')
parser.add_argument('--time_adjusted_tree', default=False, action='store_true', help='this will modify the simulated tree such that all drift lengths from root to leaf are the same')
parser.add_argument('--sadmix_tree', default=False, action='store_true', help='this will simulate trees where all admixture events are important in the sense that they expand the space of possible covariance matrices.')

#chain data collection
parser.add_argument('--summary_majority_tree', action='store_true', default=False, help='this will calculate the majority (newick) tree based on the sampled tree')
parser.add_argument('--summary_acceptance_rate', action='store_true', default=True, help='This will calculate and store summaries related to the acceptance rate')
parser.add_argument('--summary_admixture_proportion_string', action='store_true', default=True, help='this will save a string in each step indicating names and values of all admixture proportions')

#MCMCMC setup
parser.add_argument('--n', type=int, default=200, help='the number of MCMCMC flips throughout the chain.')
parser.add_argument('--m', type=int, default=50, help='the number of MCMC steps before the chain is ')
parser.add_argument('--max_temp', type=float, default=100, help='the maximum temperature used in the MCMCMC.')
parser.add_argument('--thinning_coef', type=int, default=40, help='the number of iterations between each data recording point.')
parser.add_argument('--store_permuts', action='store_true', default=False, help='If applied, the permutations from the MCMCMC flips are recorded in a file with a similar filename to the result_file')
parser.add_argument('--stop_criteria', action='store_true', default=True, help='If applied the MCMCMC will stop when the coldest chain has an effective sample size at ')
parser.add_argument('--stop_criteria_frequency', type=int, default=200000, help='This tells the frequency of checking for when the stop criteria are checked (if the stop_criteria flag is turned on)')

options=parser.parse_args()

def get_proposals(options):
    all_proposals=['deladmix', 'addadmix', 'rescale', 
                   'regraft', 'rescale_add', 'rescale_admixtures',
                   'rescale_admix_correction', 
                   'rescale_constrained', 'rescale_marginally', 
                   'sliding_regraft', 'sliding_rescale']
    all_proportions=[options.deladmix, options.addadmix, options.rescale, 
                     options.regraft, options.rescale_add, options.rescale_admix, 
                     options.rescale_admix_correction,
                     options.rescale_constrained, options.rescale_marginally, 
                     options.sliding_regraft, options.sliding_rescale]
    
    thinned_proportions=[]
    thinned_proposals=[]
    
    for proposal, proportion in zip(all_proposals, all_proportions):
        if proportion> 1e-8:
            thinned_proportions.append(proportion)
            thinned_proposals.append(proposal)
    return thinned_proportions, thinned_proposals

if options.cancel_preserve_root_distance:
    extras={'deladmix':{'preserve_root_distance':False}, 'addadmix':{'preserve_root_distance':False}}
else:
    extras={}
    
if options.no_add:
    extras['rescale_constrained']={'update_add':False}

proportions, proposals = get_proposals(options)
mp= [simple_adaptive_proposal(proposals, proportions, extras) for _ in xrange(options.MCMC_chains)]

before_added_outgroup, full_nodes, reduced_nodes=get_nodes(options.nodes, options.input_file, options.outgroup_name, options.reduce_node)
print 'before_nodes', before_added_outgroup
print 'full_nodes', full_nodes
print 'reduced_nodes', reduced_nodes

if options.prefix[-1]!='_':
    prefix=options.prefix+'_'
else:
    prefix=options.prefix
    
treemix_file=prefix+"treemix_in.txt" 
treemix_out_files=prefix+'treemix_out'


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
                          treemix_out_files=treemix_out_files,
                          ms_variance_correction=options.ms_variance_correction,
                          scale_tree_factor=options.scale_tree_factor,
                          prefix=prefix,
                          t_adjust_tree=options.time_adjusted_tree,
                          sadmix=options.sadmix_tree,
                          add_wishart_noise_to_covariance=options.wishart_noise,
                          df_of_wishart_noise_to_covariance=options.wishart_df,
                          add_file=options.starting_add)

if options.treemix_instead or options.treemix_also:
    if options.alternative_treemix_infile:
        treemix_in_file=options.alternative_treemix_infile
    else:
        treemix_in_file=treemix_file+'.gz'
    dir = os.path.dirname(__file__)
    program=os.path.join(dir,'treemixrunner.py')
    calls_to_treemix=[['python',program, 
                       '-p', str(options.treemix_processes), 
                       '-n', str(options.treemix_reps), 
                       '-i', treemix_in_file,
                       '-o', prefix+'treemix'+str(k),
                       '-m', str(k)] for k in options.treemix_no_admixtures]
    from subprocess import call
    for c in calls_to_treemix:
        call(c)
    
    if not options.treemix_also:
        import sys
        sys.exit()

no_pops=len(reduced_nodes)

if options.estimate_bootstrap_df:
    assert 6 in options.covariance_pipeline, 'Can not estimate the degrees of freedom without SNP data.'
    reduce_also= (8 in options.covariance_pipeline)
    wishart_df=estimate_degrees_of_freedom(treemix_file, 
                                           bootstrap_blocksize=options.bootstrap_blocksize, 
                                           reduce_also=reduce_also,
                                           reducer=options.reduce_node,
                                           no_bootstrap_samples=options.no_bootstrap_samples,
                                           outfile=treemix_out_files)
else:
    wishart_df=options.wishart_df
    
with open(prefix+'wishart_DF.txt', 'w') as f:
    f.write(str(wishart_df))

if not options.starting_trees:
    starting_trees=map(str, [no_pops]*options.MCMC_chains)
else:
    starting_trees=options.starting_trees

scaling=1
if 5 in options.covariance_pipeline:
    scaling*=options.scale_tree_factor    
if 9 in options.covariance_pipeline:
    scaling*=covariance[1]
    
starting_trees=get_starting_trees(starting_trees, options.MCMC_chains, options.random_start, nodes=reduced_nodes, add=options.starting_add, scaling=scaling)

summary_verbose_scheme, summaries=get_summary_scheme(majority_tree=options.summary_majority_tree, 
                                          full_tree=True, #can not think of a moment where you don't want this.
                                          proposals=mp[0], 
                                          acceptance_rate_information=options.summary_acceptance_rate,
                                          admixture_proportion_string=options.summary_admixture_proportion_string,
                                          no_chains=options.MCMC_chains)

sim_lengths=[options.m]*options.n
if 9 in options.covariance_pipeline:
    posterior, multiplier=initialize_posterior(covariance[0], M=options.wishart_df, p=options.p, 
                                               use_skewed_distr=options.sap_analysis, 
                                               multiplier=covariance[1], nodes=reduced_nodes,
                                               use_uniform_prior=options.uniform_prior)
else:
    posterior=initialize_posterior(covariance, M=options.wishart_df, p=options.p, 
                                   use_skewed_distr=options.sap_analysis, 
                                   multiplier=None, nodes=reduced_nodes,
                                   use_uniform_prior=options.uniform_prior)
    multiplier=None
    
print 'starting_trees', starting_trees
print 'posterior', posterior
print 'summaries',summaries
print 'temperature_scheme', fixed_geometrical(options.max_temp,options.MCMC_chains)
print 'summary_verbose_scheme',summary_verbose_scheme
print 'sim_lengths',sim_lengths
print 'int(options.thinning_coef)',int(options.thinning_coef)
print 'mp',mp
print 'options.MCMC_chains',options.MCMC_chains
print 'multiplier', multiplier
print 'result_file', options.result_file
print 'options.store_permuts', options.store_permuts

if options.stop_criteria:
    sc=stop_criteria(frequency=options.stop_criteria_frequency, outfile=prefix+'stop_criteria.txt')
else:
    sc=None

def f():
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
    
def g():
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
    
if options.profile:
    import cProfile
    cProfile.run('g()')
elif options.MCMC_chains==1:
    g()
else:
    res=f()
    if options.store_permuts:
        save_permuts_to_csv(res, get_permut_filename(options.result_file))


        
