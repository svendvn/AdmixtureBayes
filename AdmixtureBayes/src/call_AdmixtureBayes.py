from argparse import ArgumentParser
from meta_proposal import simple_adaptive_proposal
    
rser = ArgumentParser(usage='pipeline for Admixturebayes', version='1.0.0')
parser.add_argument('--wishart_df', type=float, default=10000.0, help='degrees of freedom to run under if bootstrap-mle of this number is declined.')
parser.add_argument('--estimate_bootstrap_df', action='store_true', default=False, help= 'if declared, the program will estimate the degrees of freedom in the wishart distribution with a bootstrap sample.')
parser.add_argument('--covariance_pipeline', nargs='+', type=int, default=[6,7,8,9], help='skewed admixture proportion prior in the simulated datasets')
parser.add_argument('--sap_analysis',  action='store_true',default=False, help='skewed admixture proportion prior in the analysis')
parser.add_argument('--input_file', type=str, default='', help='the input file of the pipeline. Its type should match the first argument of covariance_pipeline. 6= treemix file, 7-9=covariance file')
parser.add_argument('--result_file', type=str, default='result_mc3.csv', help='file to save results in')
parser.add_argument('--outgroup_name', type=str, default='default_outgroup', help='the name of the outgroup in the dataset.')

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
parser.add_argument('--random_start', action='store_true', default=False, help='If supplied, the starting trees will be simulated from the prior (unless the starting trees')

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



        
