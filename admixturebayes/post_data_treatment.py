'''
Created on May 3, 2017

@author: svendvn
'''
import summary

from tree_statistics import *
from posterior import initialize_posterior
from Rtree_to_covariance_matrix import make_covariance
import Rtree_operations
import tree_statistics
import analyse_results
from tree_plotting import plot_graph, plot_as_directed_graph


def extend_summary(filename, df=10, resfile='test.csv', other_trees=[], other_trees_tree_files=[]):
    ad=file_to_list_of_rows(filename)
    print ad
    _,names,_,values= ad
    print values[11]
    true_tree=identifier_to_tree_clean(values[11])
    plot_as_directed_graph(true_tree)
    summaries=get_summaries(true_tree, df)
    analyse_results.generate_summary_csv(summaries, reference_tree=true_tree, filename=resfile)
    
    

def file_to_list_of_rows(filename):
    tr=[]
    with open(filename, 'r') as f:
        for l in f.readlines():
            tr.append(l.rstrip().split(','))
    return zip(*tr)

        

def get_summaries(true_tree, df=10):
    m=make_covariance(true_tree)
    posterior=initialize_posterior(m, df)
    summaries=[summary.s_variable_recalculated('posterior', output='double', pks_function=posterior),
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
               summary.s_variable('proposal_type', output='string'),
               summary.s_variable('sliding_regraft_adap_param', output='double_missing'),
               summary.s_variable('rescale_adap_param', output='double_missing'),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable_recalculated(s,output='double', pks_function=posterior) for s in ['prior','branch_prior','no_admix_prior','top_prior']]
    return summaries

if __name__=='__main__':
    extend_summary('../../../Documents/9_2b/summaries.csv', df=10000, resfile='../../../Documents/9_2b/summaries_df10000.csv')
    