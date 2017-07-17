from argparse import ArgumentParser
from meta_proposal import simple_adaptive_proposal
from construct_starting_trees_choices import get_starting_trees
from construct_variance_choices import get_covariance
from construct_nodes_choices import get_nodes
    
parser = ArgumentParser(usage='pipeline for Admixturebayes', version='1.0.0')
parser.add_argument('--wishart_df', type=float, default=10000.0, help='degrees of freedom to run under if bootstrap-mle of this number is declined.')
parser.add_argument('--estimate_bootstrap_df', action='store_true', default=False, help= 'if declared, the program will estimate the degrees of freedom in the wishart distribution with a bootstrap sample.')
parser.add_argument('--covariance_pipeline', nargs='+', type=int, default=[6,7,8,9], help='skewed admixture proportion prior in the simulated datasets')
parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
parser.add_argument('--input_file', type=str, default='', help='the input file of the pipeline. Its type should match the first argument of covariance_pipeline. 6= treemix file, 7-9=covariance file')
parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in')
parser.add_argument('--outgroup_name', type=str, default='', help='the name of the outgroup that should be added to a simulated dataset.')
parser.add_argument('--reduce_node', type=str, default='', help='the name of the population that should be made root.')
parser.add_argument('--nodes', type=str, nargs='+', default=[''], help= 'list of nodes of the populations or the filename of a file where the first line contains all population names. If unspecified the first line of the input_file will be used.')

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




options=parser.parse_args()

def get_proposals(options):
    all_proposals=['deladmix', 'addadmix', 'rescale', 
                   'regraft', 'rescale_add', 'rescale_admix', 
                   'rescale_constrained', 'rescale_marginally', 
                   'sliding_regraft', 'sliding_rescale']
    all_proportions=[options.delamix, options.addadmix, options.rescale, 
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
mp= [simple_adaptive_proposal(proposals, proportions) for _ in options.MCMC_chains]

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





        
