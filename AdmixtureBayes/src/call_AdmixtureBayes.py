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


parser = ArgumentParser(usage='pipeline for Admixturebayes', version='1.0.0')
parser.add_argument('--wishart_df', type=float, default=10000.0, help='degrees of freedom to run under if bootstrap-mle of this number is declined.')
parser.add_argument('--p', type=float, default=0.5, help= 'the geometrical parameter in the prior. The formula is p**x(1-p)')
parser.add_argument('--estimate_bootstrap_df', action='store_true', default=False, help= 'if declared, the program will estimate the degrees of freedom in the wishart distribution with a bootstrap sample.')
parser.add_argument('--covariance_pipeline', nargs='+', type=int, default=[6,7,8,9], help='skewed admixture proportion prior in the simulated datasets')
parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
parser.add_argument('--input_file', type=str, default='', help='the input file of the pipeline. Its type should match the first argument of covariance_pipeline. 6= treemix file, 7-9=covariance file')
parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in')
parser.add_argument('--outgroup_name', type=str, default='', help='the name of the outgroup that should be added to a simulated dataset.')
parser.add_argument('--reduce_node', type=str, default='', help='the name of the population that should be made root.')
parser.add_argument('--nodes', type=str, nargs='+', default=[''], help= 'list of nodes of the populations or the filename of a file where the first line contains all population names. If unspecified the first line of the input_file will be used. If no input file is found, there will be used standard s1,..,sn.')

parser.add_argument('--deladmix', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--addadmix', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--regraft', type=float, default=0, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_add', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_admix', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_constrained', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--rescale_marginally', type=float, default=0, help='this states the frequency of the proposal type')
parser.add_argument('--sliding_regraft', type=float, default=1, help='this states the frequency of the proposal type')
parser.add_argument('--sliding_rescale', type=float, default=0, help='this states the frequency of the proposal type')

parser.add_argument('--MCMC_chains', type=int, default=8, help='The number of chains to run the MCMCMC with.')
parser.add_argument('--starting_trees', type=str, nargs='+', default=[], help='filenames of trees to start in. If empty the trees will either be simulated with the flag --random_start or the so-called trivial tree')
parser.add_argument('--random_start', action='store_true', default=False, help='If supplied, the starting trees will be simulated from the prior (unless the starting trees are specified)')

#tree simulation
parser.add_argument('--p_sim', type=float, default=.5, help='the parameter of the geometric distribution in the distribution to simulate the true tree from.')
parser.add_argument('--popsize', type=int, default=20, help='the number of genomes sampled from each population.')
parser.add_argument('--nreps', type=int, default=100, help='How many pieces of size 500 kb should be simualted')
parser.add_argument('--treemix_file', type=str, default='tmp.treemix_in', help= 'the filename of the intermediate step that contains the ms output.')
parser.add_argument('--ms_variance_correction', default=False, action='store_true', help= 'Should the empirical covariance matrix be adjusted for finite sample size.')
parser.add_argument('--scale_tree_factor', type=float, default=0.05, help='The scaling factor of the simulated trees to make them less vulnerable to the fixation effect.')
parser.add_argument('--skewed_admixture_prior_sim', default=False, action='store_true', help='the prior tree is simulated with an uneven prior on the admixture proportions')

#chain data collection
parser.add_argument('--summary_majority_tree', action='store_true', default=False, help='this will calculate the majority (newick) tree based on the sampled tree')
parser.add_argument('--summary_acceptance_rate', action='store_true', default=False, help='This will calculate and store summaries related to the acceptance rate')
parser.add_argument('--summary_admixture_proportion_string', action='store_true', default=False, help='this will save a string in each step indicating names and values of all admixture proportions')

#MCMCMC setup
parser.add_argument('--n', type=int, default=20000, help='the number of MCMCMC flips throughout the chain.')
parser.add_argument('--m', type=int, default=50, help='the number of MCMC steps before the chain is ')
parser.add_argument('--max_temp', type=float, default=800, help='the maximum temperature used in the MCMCMC.')
parser.add_argument('--thinning_coef', type=int, default=40, help='the number of iterations between each data recording point.')
parser.add_argument('--store_permuts', action='store_true', default=False, help='If applied, the permutations from the MCMCMC flips are recorded in a file with a similar filename to the result_file')

options=parser.parse_args()

def get_proposals(options):
    all_proposals=['deladmix', 'addadmix', 'rescale', 
                   'regraft', 'rescale_add', 'rescale_admixtures', 
                   'rescale_constrained', 'rescale_marginally', 
                   'sliding_regraft', 'sliding_rescale']
    all_proportions=[options.deladmix, options.addadmix, options.rescale, 
                     options.regraft, options.rescale_add, options.rescale_admix, 
                     options.rescale_constrained, options.rescale_marginally, 
                     options.sliding_regraft, options.sliding_rescale]
    
    thinned_proportions=[]
    thinned_proposals=[]
    
    for proposal, proportion in zip(all_proposals, all_proportions):
        if proportion> 1e-8:
            thinned_proportions.append(proportion)
            thinned_proposals.append(proposal)
    return thinned_proportions, thinned_proposals

proportions, proposals = get_proposals(options)
mp= [simple_adaptive_proposal(proposals, proportions) for _ in xrange(options.MCMC_chains)]

before_added_outgroup, full_nodes, reduced_nodes=get_nodes(options.nodes, options.input_file, options.outgroup_name, options.reduce_node)



covariance=get_covariance(options.covariance_pipeline, 
                          options.input_file, 
                          full_nodes=full_nodes, 
                          skewed_admixture_prior_sim=options.skewed_admixture_prior_sim, 
                          p=options.p_sim, 
                          outgroup_name=options.outgroup_name,
                          reduce_covariance_node=options.reduce_node,
                          sample_per_pop=options.popsize,
                          nreps=options.nreps,
                          treemix_file=options.treemix_file,
                          ms_variance_correction=options.ms_variance_correction,
                          scale_tree_factor=options.scale_tree_factor)

no_pops=len(reduced_nodes)

if not options.starting_trees:
    starting_trees=map(str, [no_pops]*options.MCMC_no_chains)
else:
    starting_trees=options.starting_trees
    
starting_trees=get_starting_trees(starting_trees, options.random_start)

summary_verbose_scheme, summaries=get_summary_scheme(majority_tree=options.summary_majority_tree, 
                                          full_tree=True, #can not think of a moment where you don't want this.
                                          proposals=mp[0], 
                                          acceptance_rate_information=options.summary_acceptance_rate,
                                          admixture_proportion_string=options.summary_admixture_proportion_string,
                                          no_chains=options.MCMC_no_chains)

sim_lengths=[options.m]*options.n
if 9 in options.covariance_pipeline:
    posterior, multiplier=initialize_posterior(covariance[0], M=options.wishart_df, p=options.m, use_skewed_distr=options.sap_analysis, multiplier=covariance[1])
else:
    posterior=initialize_posterior(covariance[0], M=options.wishart_df, p=options.m, use_skewed_distr=options.sap_analysis, multiplier=None)
    multiplier=None
    

res=MCMCMC(starting_trees=starting_trees, 
       posterior_function= posterior,
       summaries=summaries, 
       temperature_scheme=fixed_geometrical(options.max_temp,options.MCMC_no_chains), 
       printing_schemes=summary_verbose_scheme, 
       iteration_scheme=sim_lengths, 
       overall_thinnings=int(options.thinning_coef), 
       proposal_scheme= mp, 
       cores=options.MCMC_no_chains, 
       no_chains=options.MCMC_no_chains,
       multiplier=multiplier,
       result_file=options.result_file,
       store_permuts=options.store_permuts)

if options.store_permuts:
    save_permuts_to_csv(res, get_permut_filename(options.result_file))


        
