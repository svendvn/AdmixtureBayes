from argparse import ArgumentParser
from downstream_analysis_tool import (thinning, iterate_over_output_file, always_true, make_Rtree, make_full_tree, read_true_values,
                                      make_Rcovariance, cov_truecov, topology_identity,get_pops,compare_pops,extract_number_of_sadmixes, 
                                      read_one_line,summarize_all_results, create_treemix_csv_output, topology, float_mean, mode,
                                      create_treemix_sfull_tree_csv_output, subgraph, subsets, thinning_on_admixture_events)
from subgraphing import read_subgraphing_dict, get_most_likely_subgraphs_list, save_top_subgraphs
from numpy import mean
from copy import deepcopy
from collections import Counter
from Rtree_operations import get_leaf_keys


possible_summaries={'Rtree': make_Rtree,
                    'full_tree':make_full_tree,
                    'Rcov':make_Rcovariance,
                    'cov_dist':cov_truecov,
                    'topology':topology,
                    'subgraph':subgraph,
                    'subsets':subsets,
                    'top_identity':topology_identity,
                    'pops':get_pops,
                    'set_differences':compare_pops,
                    'no_sadmixes':extract_number_of_sadmixes
                    }


parser = ArgumentParser(usage='pipeline for post analysis', version='1.0.0')

parser.add_argument('--test_run', default=False, action='store_true', help='will overwrite everything and run a test function')
parser.add_argument('--input_file', default='result_mc3.csv', type=str, help='The input file that should contain a column named tree or the option no_header should be turned on in which case every line is assumed to hold just one tree')
parser.add_argument('--nodes', default='', type=str, help='file where the first line is the leaf nodes')
parser.add_argument('--subnodes', default=[], type=str, nargs='+', help='The subnodes to do the analysis on. If not declared, the analysis will be done on the full nodes.')
parser.add_argument('--no_sort', default=False, action='store_true', help='often the tree is sorted according to the leaf names. no_sort willl assumed that they are not sorted according to this but sorted according to ')
parser.add_argument('--burn_in_fraction', default=0.5, type=float, help='the proportion of the rows that are discarded as burn in period')
parser.add_argument('--total', default=886, type=int, help='an upper limit on the number of rows to reduce computational pressure')
parser.add_argument('--prefix', default='sletmig/', type=str,help='place to put the temporary files')
parser.add_argument('--result_file', default='sletmig/row.txt', type=str,help='the result file')
parser.add_argument('--outgroup_name', default='out', type=str, help='name of the outgroup')
parser.add_argument('--min_w', default=0.0, type=float, help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')
parser.add_argument('--use_cols', default=['tree','add','layer','no_admixes'], type=str, nargs='+', help='The columns to load from the input file')
parser.add_argument('--constrain_number_of_admixes', default='', type=str, choices=['','true_val']+map(str, range(21)), help='The number of admixture events that there are constrained on in the data set. If negative there are no constraints')
parser.add_argument('--constrain_number_of_effective_admixes', default='',choices=['','true_val']+map(str, range(21)), type=str, help='The number of effective_admixture events that there are constrained on in the data set. If negative there are no constraints.')
parser.add_argument('--choice_if_no_thinned_graphs', default='error', choices=['error', 'nearest_admixture_events'], help='If the thinning leaves no graphs left, this is what will be done in stead. error will throw an error and nearest_admixture_events will expand the band of allowed number of admixture events(if the chain has been thinned on number of admixture events).')
parser.add_argument('--constrain_sadmix_trees', default=False, action='store_true', help='this will remove the ')
parser.add_argument('--summaries', default=['Rtree','Rcov','cov_dist','topology','top_identity','pops','set_differences'], choices=possible_summaries.keys(),nargs='*', type=str, help='The summaries to calculate')
parser.add_argument('--save_summaries', default=['no_admixes','add','cov_dist','top_identity','set_differences'], nargs='*', type=str, help='The list of summaries to save')
parser.add_argument('--summary_summaries', default=['mean'], nargs='*', type=str, help='How each list is summarized as a single, numerical value. If it doesnt have the same length as save summaries the arguments will be repeated until it does')
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
parser.add_argument('--subgraph_file', default='', type=str, help='file where each line has a space separated list of leaf labels to calculate subtrees from. If a double underscore(__) occurs, it means that the following two arguments are max number of sub topologies and total posterior probability.')

options= parser.parse_args()

if options.subnodes:
    if options.outgroup_name in options.subnodes:
        subnodes_with_outgroup=options.subnodes
        subnodes_wo_outgroup=deepcopy(options.subnodes)
        subnodes_wo_outgroup.remove(options.outgroup_name)
    else:
        subnodes_with_outgroup=deepcopy(options.subnodes)+[options.outgroup_name]
        subnodes_wo_outgroup=options.subnodes
else:
    subnodes_with_outgroup=[]
    subnodes_wo_outgroup=[]
    

outp=read_true_values(true_scaled_tree=options.true_scaled_tree, 
                      true_tree=options.true_tree,
                      true_add=options.true_add,
                      true_covariance_reduced=options.true_covariance_reduced,
                      true_covariance_and_multiplier=options.true_covariance_and_multiplier,
                      true_no_admix=options.true_no_admix,
                      subnodes_with_outgroup=subnodes_with_outgroup,
                      subnodes_wo_outgroup=subnodes_wo_outgroup)
true_scaled_tree, true_tree, true_add, true_covariance_reduced, (true_covariance_scaled,true_multiplier), true_no_admix, _, _, _=outp
outp=read_true_values(true_covariance_reduced=options.emp_covariance_reduced,
                      true_covariance_and_multiplier=options.emp_covariance_and_multiplier,
                      true_m_scale=options.emp_m_scale,
                      subnodes_with_outgroup=subnodes_with_outgroup,
                      subnodes_wo_outgroup=subnodes_wo_outgroup)
_, _, _, emp_covariance_reduced, (emp_covariance_scaled,multiplier), _, emp_m_scale, vc, df=outp

if options.treemix_post_analysis:
    if not options.treemix_full_tree:
        outp=read_true_values(true_tree=options.treemix_tree,
                          true_add=options.treemix_add,
                      subnodes_with_outgroup=subnodes_with_outgroup,
                      subnodes_wo_outgroup=subnodes_wo_outgroup)
        _, treemix_tree, treemix_add, _, _, _, _, _, _=outp
        create_treemix_csv_output(treemix_tree,treemix_add*multiplier, emp_m_scale, options.treemix_csv_output)
    elif options.treemix_full_tree:
        outp=read_true_values(true_scaled_tree=options.treemix_full_tree,
                      subnodes_with_outgroup=subnodes_with_outgroup,
                      subnodes_wo_outgroup=subnodes_wo_outgroup)
        full_treemix_tree, _, _, _, _, _, _, _, _=outp
        create_treemix_sfull_tree_csv_output(full_treemix_tree, emp_m_scale, options.treemix_csv_output)
        full_nodes=sorted(get_leaf_keys(full_treemix_tree))
            

if options.constrain_number_of_admixes:
    if options.constrain_number_of_admixes=='true_val':
        thinner=thinning_on_admixture_events(burn_in_fraction=options.burn_in_fraction, total=options.total, no_admixes=true_no_admix, if_no_trees=options.choice_if_no_thinned_graphs)
    else:
        thinner=thinning_on_admixture_events(burn_in_fraction=options.burn_in_fraction, total=options.total, no_admixes=options.constrain_number_of_admixes,if_no_trees=options.choice_if_no_thinned_graphs)
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
possible_summary_summaries={'mean':float_mean}

print 'subnodes_wo_outgroup', subnodes_wo_outgroup

if 'Rtree' in options.summaries:
    row_sums.append(possible_summaries['Rtree'](deepcopy(nodes),options.constrain_sadmix_trees, subnodes=subnodes_wo_outgroup))
    name_to_rowsum_index('Rtree')
if 'full_tree' in options.summaries:
    row_sums.append(possible_summaries['full_tree'](add_multiplier=1.0/multiplier, 
                                                    outgroup_name=options.outgroup_name, 
                                                    remove_sadtrees=options.constrain_sadmix_trees, 
                                                    subnodes=subnodes_with_outgroup))
    name_to_rowsum_index('full_tree')
if options.subnodes:
    nodes=subnodes_wo_outgroup
    full_nodes=sorted(list(set(nodes[:]+[options.outgroup_name])))
if 'subgraph' in options.summaries:
    subgraph_dicts=read_subgraphing_dict(options.subgraph_file, types=['full'])
    for dic in subgraph_dicts:
        skeys=dic['subgraph_keys']
        identifier='.'.join(skeys)
        code='subgraph_'+identifier
        sum_func=possible_summaries['subgraph'](skeys,identifier)
        row_sums.append(sum_func)
        name_to_rowsum_index(code)
        options.save_summaries.append(code)
        options.summary_summaries.append(code)
        possible_summary_summaries[code]=sum_func.summarise
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
if 'subsets' in options.summaries:
    subgraph_dicts=read_subgraphing_dict(options.subgraph_file, types=['topological'])
    for dic in subgraph_dicts:
        skeys=dic['subgraph_keys']
        identifier='.'.join(skeys)
        code='subsets_'+identifier
        sum_func=possible_summaries['subsets'](skeys,identifier)
        row_sums.append(sum_func)
        name_to_rowsum_index(code)
        options.save_summaries.append(code)
        options.summary_summaries.append(code)
        possible_summary_summaries[code]=sum_func.summarise
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
#print 'all_results:'
#print all_results
    

def save_wrapper(filename):
    def save(listi):
        with open(filename,'w') as f:
            for ele in listi:
                f.write(str(ele))
    

    
if 'mode_topology_compare' in options.summary_summaries or 'mode_topology' in options.summary_summaries:
    def mode_topology_compare(v):
        a=mode(v)
        #print a
        return row_sums[name_to_rowsum_index['top_identity']](a)[0]['top_identity']
    def mode_topology(v):
        #print 'tops', v
        a=mode(v)
        return(a)
    possible_summary_summaries['mode_topology']=mode_topology
    possible_summary_summaries['mode_topology_compare']=mode_topology_compare
if 'mode_pops_compare' in options.summary_summaries or 'mode_pops' in options.summary_summaries:
    def mode_pops(v):
        #print 'pops',v
        v2=['-'.join(sorted(vi)) for vi in v]
        vmax_s=mode(v2)
        return vmax_s
    def mode_pops_compare(v):
        v2=['-'.join(sorted(vi)) for vi in v]
        vmax_s=mode(v2)
        vmax=vmax_s.split('-')
        #print vmax
        return row_sums[name_to_rowsum_index['set_differences']](vmax)[0]['set_differences']
    possible_summary_summaries['mode_pops_compare']=mode_pops_compare
    possible_summary_summaries['mode_pops']=mode_pops
if 'node_count' in options.summary_summaries:
    def count_measure(v):
        ad=Counter([a for k in v for a in k])
        l=ad.most_common(4000)
        total=len(v)
        with open('node_counts.txt','w') as f:
            for key,count_num in l:
                f.write(key+' '+str(float(count_num)/total)+'\n')
    possible_summary_summaries['node_count']=count_measure
if 'top_pops' in options.summary_summaries:
    def write_top_pops(v):
        v2=['_'.join(sorted(vi)) for vi in v]
        ad=Counter(v2)
        l=ad.most_common(10)
        total=len(v)
        with open('top_pops.txt','w') as f:
            for n,(key,count_num) in enumerate(l):
                f.write(str(n+1)+','+str(float(count_num)/total)+','+key+'\n')
    possible_summary_summaries['top_pops']=write_top_pops
if 'subgraph' in options.summary_summaries:
    subgraph_dicts=read_subgraphing_dict(options.subgraph_file)
    def subgraphing(trees):
        for dic in subgraph_dicts:
            sgraphs=get_most_likely_subgraphs_list(strees=trees, nodes=nodes, subgraph_keys=dic['subgraph_keys'])
            save_top_subgraphs(topologies=sgraphs, nodes=dic['subgraph_keys'], **dic)
    possible_summary_summaries['subgraph']=subgraphing
    

    
    
    


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
    

        
