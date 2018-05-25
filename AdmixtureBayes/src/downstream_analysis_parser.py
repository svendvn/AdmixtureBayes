from argparse import ArgumentParser
from downstream_analysis_tool import thinning, iterate_over_output_file, always_true, make_Rtree, make_full_tree, read_true_values


possible_summaries={'Rtree': make_Rtree,
                    'full_tree':make_full_tree,
                    'Rcov':make_Rcovariance,
                    'cov_dist':cov_truecov,
                    'top_identity':topology_identity,
                    'pops':get_pops,
                    'set_differences':compare_pops,
                    'no_sadmixes':extract_number_of_admixes
                    }

parser = ArgumentParser(usage='pipeline for post analysis', version='1.0.0')

parser.add_argument('--test_run', default=False, action='store_true', help='will overwrite everything and run a test function')
parser.add_argument('--input_file', default='result_mc3.csv', type=str, help='The input file that should contain a column named tree or the option no_header should be turned on in which case every line is assumed to hold just one tree')
parser.add_argument('--nodes', default='', type=str, help='file where the first line is the leaf nodes')
parser.add_argument('--no_sort', default=False, action='store_true', help='often the tree is sorted according to the leaf names. no_sort willl assumed that they are not sorted according to this but sorted according to ')
parser.add_argument('--burn_in_fraction', default=0.0, type=float, help='the proportion of the rows that are discarded as burn in period')
parser.add_argument('--total', default=500, type=int, help='an upper limit on the number of rows to reduce computational pressure')
parser.add_argument('--prefix', default='sletmig/', type=str,help='place to put the temporary files')
parser.add_argument('--result_file', default='sletmig/row.txt', type=str,help='the result file')
parser.add_argument('--outgroup_name', default='out', type=str, help='name of the outgroup')
parser.add_argument('--min_w', default=0.0, type=float, help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')
parser.add_argument('--use_cols', default=['tree','add','layer','no_admixes'], type=str, nargs='+', help='The columns to load from the input file')
parser.add_argument('--constrain_number_of_admixes', default='', type=str, choices=['','true_val']+map(str, range(21)), help='The number of admixture events that there are constrained on in the data set. If negative there are no constraints')
parser.add_argument('--constrain_number_of_effective_admixes', default='',choices=['','true_val']+map(str, range(21)), type=str, help='The number of effective_admixture events that there are constrained on in the data set. If negative there are no constraints.')
parser.add_argument('--summaries', default=['Rtree','full_tree','Rcov','cov_dist'], choices=possible_summaries.keys(),nargs='+', type=str, help='The summaries to calculate')
parser.add_argument('--save_summaries', default=['no_admixes'], nargs='+', type='str', help='The list of summaries to save')
parser.add_argument('--summary_summaries', default=['mean'], nargs='+', type='str', help='How each list is summarized as a single, numerical value. If it doesnt have the same length as save summaries the arguments will be repeated until it does')
parser.add_argument('--true_scaled_tree',  type=str, default='')
parser.add_argument('--true_tree',  type=str, default='')
parser.add_argument('--true_add',  type=str, default='')
parser.add_argument('--true_covariance_reduced',  type=str, default='')
parser.add_argument('--true_covariance_and_multiplier',  type=str, default='')
parser.add_argument('--true_no_admix',  type=str, default='')
parser.add_argument('--true_m_scale',  type=str, default='')
parser.add_argument('--variance_correction',  type=str, default='')
parser.add_argument('--df',  type=str, default='')

options= parser.parse_args()

outp=read_true_values(true_scaled_tree=options.true_scaled_tree, 
                      true_tree=options.true_tree,
                      true_add=options.true_add,
                      true_covariance_reduced=options.true_covariance_reduced,
                      true_covariance_and_multiplier=options.true_covariance_and_multiplier,
                      true_no_admix=options.true_no_admix,
                      true_m_scale=options.true_m_scale,
                      variance_correction=options.variance_correction,
                      df=options.df)
true_scaled_tree, true_tree, true_add, true_covariance_reduced, (true_covariance_scaled,multiplier), true_no_admix, true_m_scale, vc, df=outp

if options.constrain_number_of_admixes:
    if options.constrain_number_of_admixes=='true_val':
        thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.total, no_admixes=true_no_admix)
    else:
        thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.total, no_admixes=options.constrain_number_of_admixes)
else:
    thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.total)

nodes=read_one_line(options.nodes)
if not options.no_sort:
    nodes=sorted(nodes)

row_sums=[]


    

    
all_results=iterate_over_output_file(outfile, 
                                     cols=options.use_cols, 
                                     pre_thin_data_set_function=thinner, 
                                     while_thin_data_set_functions=always_true,
                                     row_summarize_functions=[],
                                     thinned_d_dic=[],
                                     full_summarize_functions=[],
                                     **constant_kwargs)

