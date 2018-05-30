from argparse import ArgumentParser
from downstream_analysis_tool import (thinning, iterate_over_output_file, always_true, make_Rtree, make_full_tree, read_true_values,
                                      make_Rcovariance, cov_truecov, topology_identity,get_pops,compare_pops,extract_number_of_sadmixes, 
                                      read_one_line,summarize_all_results, create_treemix_csv_output, topology, float_mean, mode,
                                      create_treemix_sfull_tree_csv_output)
from numpy import mean
from copy import deepcopy
from collections import Counter
from Rtree_operations import get_leaf_keys


possible_summaries={'Rtree': make_Rtree,
                    'full_tree':make_full_tree,
                    'Rcov':make_Rcovariance,
                    'cov_dist':cov_truecov,
                    'topology':topology,
                    'top_identity':topology_identity,
                    'pops':get_pops,
                    'set_differences':compare_pops,
                    'no_sadmixes':extract_number_of_sadmixes
                    }


parser = ArgumentParser(usage='pipeline for post analysis', version='1.0.0')

parser.add_argument('--test_run', default=False, action='store_true', help='will overwrite everything and run a test function')
parser.add_argument('--input_file', default='result_mc3.csv', type=str, help='The input file that should contain a column named tree or the option no_header should be turned on in which case every line is assumed to hold just one tree')
parser.add_argument('--nodes', default='', type=str, help='file where the first line is the leaf nodes')
parser.add_argument('--no_sort', default=False, action='store_true', help='often the tree is sorted according to the leaf names. no_sort willl assumed that they are not sorted according to this but sorted according to ')
parser.add_argument('--burn_in_fraction', default=0.2, type=float, help='the proportion of the rows that are discarded as burn in period')
parser.add_argument('--total', default=886, type=int, help='an upper limit on the number of rows to reduce computational pressure')
parser.add_argument('--prefix', default='sletmig/', type=str,help='place to put the temporary files')
parser.add_argument('--result_file', default='sletmig/row.txt', type=str,help='the result file')
parser.add_argument('--outgroup_name', default='out', type=str, help='name of the outgroup')
parser.add_argument('--min_w', default=0.0, type=float, help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')
parser.add_argument('--use_cols', default=['tree','add','layer','no_admixes'], type=str, nargs='+', help='The columns to load from the input file')
parser.add_argument('--constrain_number_of_admixes', default='', type=str, choices=['','true_val']+map(str, range(21)), help='The number of admixture events that there are constrained on in the data set. If negative there are no constraints')
parser.add_argument('--constrain_number_of_effective_admixes', default='',choices=['','true_val']+map(str, range(21)), type=str, help='The number of effective_admixture events that there are constrained on in the data set. If negative there are no constraints.')
parser.add_argument('--constrain_sadmix_trees', default=False, action='store_true', help='this will remove the ')
parser.add_argument('--summaries', default=['Rtree','Rcov','cov_dist','topology','top_identity','pops','set_differences'], choices=possible_summaries.keys(),nargs='+', type=str, help='The summaries to calculate')
parser.add_argument('--save_summaries', default=['no_admixes','add','cov_dist','top_identity','set_differences'], nargs='+', type=str, help='The list of summaries to save')
parser.add_argument('--summary_summaries', default=['mean'], nargs='+', type=str, help='How each list is summarized as a single, numerical value. If it doesnt have the same length as save summaries the arguments will be repeated until it does')
parser.add_argument('--true_scaled_tree',  type=str, default='')
parser.add_argument('--true_tree',  type=str, default='')
parser.add_argument('--true_add',  type=str, default='')
parser.add_argument('--true_covariance_reduced',  type=str, default='')
parser.add_argument('--true_covariance_and_multiplier',  type=str, default='')
parser.add_argument('--true_no_admix',  type=str, default='')
parser.add_argument('--emp_m_scale',  type=str, default='')
parser.add_argument('--emp_variance_correction',  type=str, default='')
parser.add_argument('--emp_df',  type=str, default='')
parser.add_argument('--emp_covariance_and_multiplier', default='', type=str)
parser.add_argument('--emp_covariance_reduced', default='', type=str)
parser.add_argument('--treemix_post_analysis', action='store_true', default=False, help='this will convert the treemix input fil ../../../../Dropbox/Bioinformatik/AdmixtureBayes/test_final_grid/ai_2_5true/_true_tree.txtes into a suitable csv file for ')
parser.add_argument('--treemix_tree', default='', type=str, help='')
parser.add_argument('--treemix_add', default='', type=str, help='')
parser.add_argument('--treemix_full_tree', default='')
parser.add_argument('--treemix_csv_output', default='treemix.csv', type=str, help='')

options= parser.parse_args()

outp=read_true_values(true_scaled_tree=options.true_scaled_tree, 
                      true_tree=options.true_tree,
                      true_add=options.true_add,
                      true_covariance_reduced=options.true_covariance_reduced,
                      true_covariance_and_multiplier=options.true_covariance_and_multiplier,
                      true_no_admix=options.true_no_admix)
true_scaled_tree, true_tree, true_add, true_covariance_reduced, (true_covariance_scaled,true_multiplier), true_no_admix, _, _, _=outp
outp=read_true_values(true_covariance_reduced=options.emp_covariance_reduced,
                      true_covariance_and_multiplier=options.emp_covariance_and_multiplier,
                      true_m_scale=options.emp_m_scale)
_, _, _, emp_covariance_reduced, (emp_covariance_scaled,multiplier), _, emp_m_scale, vc, df=outp

if options.treemix_post_analysis:
    if not options.treemix_full_tree:
        outp=read_true_values(true_tree=options.treemix_tree,
                          true_add=options.treemix_add)
        _, treemix_tree, treemix_add, _, _, _, _, _, _=outp
        create_treemix_csv_output(treemix_tree,treemix_add*multiplier, emp_m_scale, options.treemix_csv_output)
    elif options.treemix_full_tree:
        outp=read_true_values(true_scaled_tree=options.treemix_full_tree)
        full_treemix_tree, _, _, _, _, _, _, _, _=outp
        create_treemix_sfull_tree_csv_output(full_treemix_tree, emp_m_scale, options.treemix_csv_output)
        full_nodes=sorted(get_leaf_keys(full_treemix_tree))
            

if options.constrain_number_of_admixes:
    if options.constrain_number_of_admixes=='true_val':
        thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.total, no_admixes=int(true_no_admix))
    else:
        thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.total, no_admixes=int(options.constrain_number_of_admixes))
else:
    thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.total)

nodes=read_one_line(options.nodes).split()
if not options.no_sort:
    nodes=sorted(nodes)

row_sums=[]

class pointers(object):
    
    def __init__(self):
        self.count=0
        self.dic={}
        
    def __call__(self, name):
        self.dic[name]=self.count
        self.count+=1
    
    def __getitem__(self, key):
        return self.dic[key]
    
name_to_rowsum_index=pointers()
    

if 'Rtree' in options.summaries:
    row_sums.append(possible_summaries['Rtree'](deepcopy(nodes),options.constrain_sadmix_trees))
    name_to_rowsum_index('Rtree')
if 'full_tree' in options.summaries:
    row_sums.append(possible_summaries['full_tree'](add_multiplier=1.0/multiplier, outgroup_name=options.outgroup_name, remove_sadtrees=options.constrain_sadmix_trees))
    name_to_rowsum_index('full_tree')
if 'Rcov' in options.summaries:
    row_sums.append(possible_summaries['Rcov'](deepcopy(nodes), add_multiplier=1.0/multiplier))
    name_to_rowsum_index('Rcov')
if 'cov_dist' in options.summaries:
    row_sums.append(possible_summaries['cov_dist'](true_covariance_reduced))
    name_to_rowsum_index('cov_dist')
if 'topology' in options.summaries:
    row_sums.append(possible_summaries['topology'](nodes=nodes))
    name_to_rowsum_index('topology')
if 'top_identity' in options.summaries:
    row_sums.append(possible_summaries['top_identity'](true_tree, nodes=nodes))
    name_to_rowsum_index('top_identity')
if 'pops' in options.summaries:
    row_sums.append(possible_summaries['pops'](min_w=options.min_w, keys_to_include=nodes))
    name_to_rowsum_index('pops')
if 'set_differences' in options.summaries:
    row_sums.append(possible_summaries['set_differences'](true_tree, min_w=options.min_w, keys_to_include=nodes))
    name_to_rowsum_index('set_differences')
if 'no_sadmixes' in options.summaries:
    if options.constrain_number_of_effective_admixes:
        no_effective_admixes=int(options.constrain_number_of_effective_admixes)
    else:
        no_effective_admixes=None
    row_sums.append(possible_summaries['no_sadmixes'](no_effective_admixes))
    name_to_rowsum_index('no_sadmixes')
    
def save_thin_columns(d_dic):
    return {summ:d_dic[summ] for summ in options.save_summaries}
    
    
if options.treemix_post_analysis:
    if options.treemix_full_tree:
        constant_kwargs={'full_nodes':full_nodes}
    else:
        constant_kwargs={}
    all_results,_=iterate_over_output_file(options.treemix_csv_output, 
                                         cols=options.use_cols, 
                                         pre_thin_data_set_function=thinner, 
                                         while_thin_data_set_function=always_true,
                                         row_summarize_functions=row_sums,
                                         thinned_d_dic=save_thin_columns,
                                         full_summarize_functions=[],
                                         **constant_kwargs)
else:
    all_results,_=iterate_over_output_file(options.input_file, 
                                         cols=options.use_cols, 
                                         pre_thin_data_set_function=thinner, 
                                         while_thin_data_set_function=always_true,
                                         row_summarize_functions=row_sums,
                                         thinned_d_dic=save_thin_columns,
                                         full_summarize_functions=[])
    
possible_summary_summaries={'mean':float_mean}
if 'mode_topology' in options.summary_summaries:
    def mode_topology(v):
        a=mode(v)
        print a
        return row_sums[name_to_rowsum_index['top_identity']](a)[0]['top_identity']
    possible_summary_summaries['mode_topology']=mode_topology
if 'mode_pops' in options.summary_summaries:
    def mode_pops(v):
        v2=['-'.join(sorted(vi)) for vi in v]
        vmax_s=mode(v2)
        vmax=vmax_s.split('-')
        print vmax
        return row_sums[name_to_rowsum_index['set_differences']](vmax)[0]['set_differences']
    possible_summary_summaries['mode_pops']=mode_pops

    
    
    


n=len(options.save_summaries)
summary_summaries=options.summary_summaries
while len(summary_summaries)<n:#repeat arguments until the number of arguments is correct
    summary_summaries+=options.summary_summaries

summary_summaries_functions=[possible_summary_summaries[summ] for summ in summary_summaries]


summ_results=summarize_all_results(all_results, options.save_summaries, summary_summaries_functions)
res=[]
header=[]
with open(options.result_file, 'w') as f:
    for n,(summ_func_name, summ_name) in enumerate(zip(summary_summaries, options.save_summaries)):
        res.append(summ_results[n])
        header.append(summ_name+'_'+summ_func_name)
    f.write(','.join(['input_file']+header)+'\n')
    f.write(','.join([options.input_file]+map(str,res)))
    

        
