

from reduce_covariance import reduce_covariance
from Rtree_operations import get_number_of_leaves, add_outgroup
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from Rtree_to_covariance_matrix import make_covariance
import numpy as np
from generate_sadmix_trees import effective_number_of_admixes, admixes_are_sadmixes
#def effective_number_of_admixes(t):
#   return np.random.choice(3)
from calculate_covariance_distances import open_cov_file_admb
from construct_covariance_choices import read_tree
import pandas as pd
import os

from tree_statistics import admixture_sorted_unique_identifier
from Rtree_operations import get_leaf_keys
print 'imported software'
l,outfile=sys.argv[1:3]

def get_list_of_turned_topologies(trees, true_tree):
    nodes=get_leaf_keys(true_tree)
    return [admixture_sorted_unique_identifier(tree, nodes) for tree in trees], admixture_sorted_unique_identifier(true_tree, nodes)


def get_covariance(outfile):
    cov, _, mult= open_cov_file_admb(outfile, None)
    return cov, mult

def always_true(*args):
    return True

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             pre_thin_data_set_function=identity, 
                             while_thin_data_set_functions=always_true(),
                             row_summarize_functions=[],
                             thinned_d_dic=[],
                             full_summarize_functions=[],
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols)
    df= pre_thin_data_set_function(df)
    full_summs=[full_summarize_function(df) for full_summarize_function in full_summarize_functions]
    
    all_results=[]
    
    n=len(df)
    for i,r in df.iterrows():
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)
        if not while_thin_data_set_function(**d_dic):
            continue
        for row_summarize_function in row_summarize_functions:
            d_dic.update(row_summarize_function(**d_dic))
        all_results.append(thinned_d_dic(d_dic))
    return all_results


class thinning(object):
    
    def __init__(self, burn_in_fraction=None, total=None, **values_to_filter_by):
        self.burn_in_fraction=burn_in_fraction
        self.total=total
        self.values_to_filter_by=values_to_filter_by
        
    def __call__(self, df):
        #first_removing burn-in
        n=len(df)
        print 'removing burn in from dataframe with ', n, 'rows.'
        if self.burn_in_fraction is not None:
            
            df=df[int(n/2):]
        print 'burn_in operation finished. Now ', n, 'rows.'
        for column,value in self.values_to_filter_by.items():
            print 'filtering on ', column,'==',value
            df=df.loc[df[column]==value,:]
        print 'after filtering operations there are now', len(df),'rows'
        if self.total is not None:
            n=len(df)
            stepsize=max(n//total,1)
            df=df[::stepsize]
            print 'thinning complete. Now', len(df), 'rows'
        return df
    
    
def thin_dataset(total, burn_in_fraction, layer,total):
    def thinner(df):
        df
        
            
        
    
    
    
def get_posterior_trees(outfile, add_multiplier=1, nodes=None, outgroup='out', total=500, constrained=None, sadmixes_only=False):
    a=pd.read_csv(outfile, usecols=['posterior','no_admixes','tree','add'])
    if constrained is not None:
        a=a.loc[a.no_admixes==constrained,:]
    #print a
    halfway=a.shape[0]//2
    stepsize=max(halfway//total,1)
    b=a[halfway::stepsize]
    matrices=[]
    trees=[]
    no_admixes=[]
    for stree, add in zip(b['tree'], b['add']):
        #print stree
        no_nodes=len(stree.split('-')[0].split('.'))
        full_nodes=sorted(['s'+str(i) for i in range(1, no_nodes+1)])
        tree=identifier_to_tree_clean(stree, leaves=generate_predefined_list_string(full_nodes))
        if sadmixes_only and (not admixes_are_sadmixes(tree)):
            continue
        #print pretty_string(tree)
        tree= add_outgroup(tree,  inner_node_name='new_node', to_new_root_length=float(add)*add_multiplier, to_outgroup_length=0, outgroup_name=outgroup)
        no_admixes.append(effective_number_of_admixes(tree))
        trees.append(tree)
        cov=make_covariance(tree, node_keys=[outgroup]+['s'+str(i) for i in range(1,get_number_of_leaves(tree)-1+1)])
        #print cov
        matrices.append(reduce_covariance(cov, 0))
    return matrices, trees, no_admixes,halfway*2
    
def get_true_tree(outfile, nodes):
    return read_tree(outfile, nodes)
    
def results_file_to_cov_file(r):
    base='_'.join(r.split(os.sep)[-1].split('_')[:-2])
    return os.path.join(base,'_covariance_and_multiplier.txt')
def results_file_to_tree_file(r):
    base='_'.join(r.split(os.sep)[-1].split('_')[:-2])
    return os.path.join(base,'_scaled_true_tree.txt')
def extract_info(l):
    pieces=l.split(os.sep)[-1].split('.')[0].split('_')[1:]
    if pieces[-1]=='perfect' or pieces[-1]=='random':
        pieces=pieces[:-1]
    return map(int, pieces)

def extract_info(l):
    pieces=l.split(os.sep)[-1].split('.')[0].split('_')[1:]
    if pieces[-1]=='perfect' or pieces[-1]=='random':
        pieces=pieces[:-1]
    return map(int, pieces)

def get_df(fil):
       with open(fil, 'r') as f:
              return(float(f.readline()))