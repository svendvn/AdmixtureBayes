from reduce_covariance import reduce_covariance
from Rtree_operations import get_number_of_leaves, add_outgroup, get_number_of_admixes, scale_tree, scale_tree_copy
from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, 
                             identifier_file_to_tree_clean, unique_identifier_and_branch_lengths)
from Rtree_to_covariance_matrix import make_covariance, get_populations
import numpy as np
from generate_sadmix_trees import effective_number_of_admixes, admixes_are_sadmixes
#def effective_number_of_admixes(t):
#   return np.random.choice(3)
from calculate_covariance_distances import open_cov_file_admb
from construct_covariance_choices import read_tree, read_one_line
import pandas as pd
import os
from copy import deepcopy
from tree_to_data import file_to_emp_cov

from tree_statistics import admixture_sorted_unique_identifier
from Rtree_operations import get_leaf_keys
print 'imported software'

from argparse import ArgumentParser
from collections import Counter


def get_list_of_turned_topologies(trees, true_tree):
    nodes=get_leaf_keys(true_tree)
    return [admixture_sorted_unique_identifier(tree, nodes) for tree in trees], admixture_sorted_unique_identifier(true_tree, nodes)


def get_covariance(outfile):
    cov, _, mult= open_cov_file_admb(outfile, None)
    return cov, mult

def always_true(**args):
    return True

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             pre_thin_data_set_function=identity, 
                             while_thin_data_set_function=always_true,
                             row_summarize_functions=[],
                             thinned_d_dic=identity,
                             full_summarize_functions=[],
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols)
    df = df[cols]
    df= pre_thin_data_set_function(df)
    full_summs=[full_summarize_function(df) for full_summarize_function in full_summarize_functions]
    
    all_results=[]
    
    n=len(df)
    for i,r in df.iterrows():
        cont=False
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)
        if not while_thin_data_set_function(**d_dic):
            continue
        for row_summarize_function in row_summarize_functions:
            add_dic, skip=row_summarize_function(**d_dic)
            if skip:
                cont=True
                break
            d_dic.update(add_dic)
        if cont:
            continue
        all_results.append(thinned_d_dic(d_dic))
    return all_results, full_summs

def create_treemix_csv_output(tree,add,m_scale, outfile):
    if m_scale is not None:
        tree=scale_tree(tree,1.0/m_scale)
        add=add/m_scale
    with open(outfile, 'w') as f:
        f.write('tree,add'+'\n')
        f.write(unique_identifier_and_branch_lengths(tree)+','+str(add))

class make_Rtree(object):
    
    def __init__(self, nodes_to_be_sorted, remove_sadtrees=False):
        self.nodes=sorted(nodes_to_be_sorted)
        self.remove_sadtrees=remove_sadtrees
        
    def __call__(self, tree, **not_needed):
        #print tree
        #print not_needed
        Rtree=identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        if self.remove_sadtrees and (not admixes_are_sadmixes(Rtree)):
            return {'Rtree':Rtree}, True
        return {'Rtree':Rtree}, False
    
class make_full_tree(object):
    
    def __init__(self, add_multiplier=1, outgroup_name='out'):
        self.add_multiplier=add_multiplier
        self.outgroup_name=outgroup_name
        
    def __call__(self, Rtree, add, **not_needed):
        full_tree= add_outgroup(Rtree,  
                                inner_node_name='new_node', 
                                to_new_root_length=float(add)*self.add_multiplier, 
                                to_outgroup_length=0, 
                                outgroup_name=self.outgroup_name)
        return {'full_tree':full_tree}, False
    
class make_Rcovariance(object):
    
    def __init__(self, nodes, add_multiplier=1):
        self.nodes=nodes
        self.add_multiplier=add_multiplier
        
    def __call__(self, Rtree, add, **not_needed):
        Rcov=make_covariance(Rtree, node_keys=self.nodes)+float(add)*self.add_multiplier
        
        return {'Rcov':Rcov}, False
    
class cov_truecov(object):
    
    def __init__(self, true_covariance):
        self.true_covariance=true_covariance
        
    def __call__(self, Rcov, **not_needed):
        dist=np.linalg.norm(Rcov-self.true_covariance)
        #print 'Rcov', Rcov
        #print 'true_cov', self.true_covariance
        return {'cov_dist':dist}, False
    
class topology(object):
    
    def __init__(self, nodes):
        self.nodes=nodes
        
    def __call__(self, Rtree, **not_needed):
        top=admixture_sorted_unique_identifier(Rtree, leaf_order=self.nodes, not_opposite=True)
        return {'topology':top}, False
    
class topology_identity(object):
    
    def __init__(self, true_Rtree, nodes):
        self.full_scaled_turned_topology=admixture_sorted_unique_identifier(true_Rtree, leaf_order=nodes, not_opposite=True)
        
    def __call__(self, topology, **not_needed):
        ident_top=(topology==self.full_scaled_turned_topology)
        return {'top_identity':ident_top}, False
    
class get_pops(object):
    
    def __init__(self, min_w=0.0):
        self.min_w=min_w
            
    def __call__(self, Rtree, **not_needed):
        pops=get_populations(Rtree, self.min_w)
        return {'pops':pops}, False 
    
class compare_pops(object):
    def __init__(self, true_Rtree, min_w=0.0):
        self.true_pops=set(get_populations(true_Rtree, min_w))
        
    def __call__(self, pops, **not_needed):
        diffs=set(pops).symmetric_difference(self.true_pops)
        return {'set_differences':len(diffs)}, False

class set_identity(object):
    def __init__(self):
        pass
    
    def __call__(self, set_differences,**not_needed):
        ident=(set_differences==0)
        return {'set_identity':ident}
    
    
class extract_number_of_sadmixes(object):
    
    def __init__(self, filter_on_sadmixes=None):
        self.filter_on_sadmixes=filter_on_sadmixes
    
    def __call__(self, Rtree, **not_needed):
        no_sadmixes=effective_number_of_admixes(Rtree)
        if self.filter_on_sadmixes is not None and no_sadmixes!= self.filter_on_sadmixes:
            return {}, True
        return {'no_sadmixes':no_sadmixes}, False


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
            df=df[int(n*self.burn_in_fraction):]
        print 'burn_in operation finished. Now ', len(df), 'rows.'
        for column,value in self.values_to_filter_by.items():
            print 'filtering on ', column,'==',value
            df=df.loc[df[column]==value,:]
        print 'after filtering operations there are now', len(df),'rows'
        if self.total is not None:
            n=len(df)
            stepsize=max(n//self.total,1)
            df=df[::stepsize]
            print 'thinning complete. Now', len(df), 'rows'
        return df
    
def summarize_all_results(all_results, summ_names, summ_funcs):
    lists={summ_name:[] for summ_name in summ_names}
    for row_dic in all_results:
        for summ_name, summ_val in row_dic.items():
            lists[summ_name].append(summ_val)  
    res=[]
    for summ_name, summ_func in zip(summ_names, summ_funcs):
        res.append(summ_func(lists[summ_name]))
    return res
    

def read_tree(filename):
    with open(filename, 'r') as f:
        nodes=f.readline().strip().split()
        stree=f.readline().strip()
        
        
        
def read_covariance(filename, nodes):
    cov=[]
    with open(filename, 'r') as f:
        nodes=f.readline().strip().split()
        count=0
        for lin in f.readlines():
            if count==len(nodes):
                break
            cov.append(lin.split()[1:])
            count+=1
        if len(lin)>4:
            multiplier=float(lin.split('=')[1])
        else:
            multiplier=None
    
def read_true_values(true_scaled_tree='', 
                      true_tree='',
                      true_add='',
                      true_covariance_reduced='',
                      true_covariance_and_multiplier='',
                      true_no_admix='',
                      true_m_scale='',
                      true_variance_correction=None,
                      true_df=''):
    scaled_tree,tree,add,covariance_reduced,(Rcovariance,multiplier), no_admix,m_scale,vc,df=None,None,None,None,(None,None),None,None,None,None
    if true_scaled_tree:
        scaled_tree=identifier_file_to_tree_clean(true_scaled_tree)
    if true_tree:
        tree=identifier_file_to_tree_clean(true_tree)
    if true_add:
        add=float(read_one_line(true_add))
    if true_covariance_reduced:
        covariance_reduced=file_to_emp_cov(true_covariance_reduced, sort_nodes_alphabetically=True)
    if true_covariance_and_multiplier:
        Rcovariance,multiplier,vc=file_to_emp_cov(true_covariance_and_multiplier, sort_nodes_alphabetically=True, vc=true_variance_correction)
    if true_no_admix:
        no_admix=int(true_no_admix)
    if true_m_scale:
        m_scale=float(read_one_line(true_m_scale))
    if true_df:
        df=float(read_one_line(true_df))
    return scaled_tree,tree,add,covariance_reduced,(Rcovariance,multiplier), no_admix,m_scale,vc,df
            
    
def float_mean(v):
    return np.mean(map(float,v))
def mode(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]
    
    
          
