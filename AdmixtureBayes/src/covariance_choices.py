from tree_statistics import identifier_to_tree_clean
from tree_to_data import file_to_emp_cov, reduce_covariance, ms_to_treemix3, call_ms_string, tree_to_ms_command
from generate_prior_trees import simulate_number_of_admixture_events, generate_phylogeny
from Rtree_operations import add_outgroup, get_number_of_leaves
from scipy.stats import expon
from Rtree_to_covariance_matrix import make_covariance
from load_data import read_data
from copy import deepcopy





def simulate_tree_wrapper(nk_tuple, **kwargs):
    return generate_phylogeny(size= nk_tuple[0], 
                              admixes=nk_tuple[1], 
                              leaf_nodes= kwargs['reduced_nodes'], 
                              skewed_admixture_prior=kwargs['skewed_admixture_prior_sim'])
    
def add_outgroup_wrapper(tree_without_outgroup, **kwargs):
    new_length_1, new_length_2= expon.rvs(), expon.rvs()
    tree_with_outgroup= add_outgroup(tree_without_outgroup,  inner_node_name='new_node', to_new_root_length=new_length_1, to_outgroup_length=new_length_2, outgroup_name=kwargs['outgroup_name'])
    return tree_with_outgroup

def simulate_number_of_admixture_events_wrapper(n, **kwargs):
    nk_tuple=(n, simulate_number_of_admixture_events(kwargs['p']))
    return nk_tuple

def theoretical_covariance_wrapper(tree, **kwargs):
    covariance= make_covariance(tree, node_keys= kwargs['full_nodes'])
    if kwargs['add_wishart_noise_to_covariance']:
        covariance=add_wishart_noise(covariance, kwargs['df_of_wishart_noist_to_covariance'])
    return covariance
        
def add_wishart_noise(matrix, df):
    r=matrix.shape[0]
    m=wishart.rvs(df=r*df-1, scale=matrix/(r*df))
    return m

def reduce_covariance_wrapper(full_covariance,**kwargs):
    outgroup_node_index=next((i for i,s in enumerate(kwargs['full_nodes']) if s==kwargs['outgroup_name']))
    return reduce_covariance(full_covariance, outgroup_node_index)
    
def ms_simulate_wrapper(tree, **kwargs):
    no_pops= get_number_of_leaves(tree)
    ms_command=tree_to_ms_command(tree,
                       sample_per_pop=kwargs['sample_per_pop'],
                       nreps=kwargs['nreps'],
                       theta=kwargs['theta'],
                       sites=kwargs['sites'],
                       recomb_rate=kwargs['recomb_rate'])
    print ms_command
    call_ms_string(ms_command, kwargs['ms_file'])
    filename_gz=ms_to_treemix3(kwargs['ms_file'], 
                   samples_per_pop=kwargs['sample_per_pop'], 
                   no_pops=no_pops, 
                   n_reps=kwargs['nreps'], 
                   filename2=kwargs['treemix_file'])
    return filename_gz

def empirical_covariance_wrapper(snp_data_file, **kwargs):
    return read_data(snp_data_file, 
                     blocksize=kwargs['blocksize_empirical_covariance'], 
                     nodes=kwargs['full_nodes'], 
                     noss=kwargs['ms_variance_correction'])



dictionary_of_transformations={
    (1,2):simulate_number_of_admixture_events_wrapper,
    (2,3):simulate_tree_wrapper,
    (3,4):add_outgroup_wrapper,
    (3,6):theoretical_covariance_wrapper,
    (4,6):theoretical_covariance_wrapper,
    (6,7):reduce_covariance_wrapper,
    (3,5):ms_simulate_wrapper,
    (4,5):ms_simulate_wrapper,
    (5,6):empirical_covariance_wrapper
    }


def get_covariance(stages_to_go_through, input, full_nodes=None,
                   skewed_admixture_prior_sim=False,
                   p=0.5,
                   outgroup_name='outgroup_name',
                   add_wishart_noise_to_covariance=False,
                   df_of_wishart_noise_to_covariance=1000,
                   reduce_covariance_node=None,
                   sample_per_pop=50, nreps=2, 
                   theta=0.4, sites=500000, recomb_rate=1,
                   ms_file='tmp.txt',
                   treemix_file='tmp.treemix_in',
                   blocksize_empirical_covariance=100,
                   ms_variance_correction=False):
    
    kwargs={}
    kwargs['skewed_admixture_prior_sim']=skewed_admixture_prior_sim
    kwargs['p']=p
    kwargs['outgroup_name']=outgroup_name
    kwargs['add_wishart_noise_to_covariance']=add_wishart_noise_to_covariance
    kwargs['df_of_wishart_noise_to_covariance']=df_of_wishart_noise_to_covariance
    kwargs['full_nodes']=full_nodes
    if reduce_covariance_node is None:
        reduce_covariance_node=outgroup_name
    kwargs['reduce_covariance_node']=reduce_covariance_node
    reduced_nodes=deepcopy(full_nodes)
    if outgroup_name in reduced_nodes:
        reduced_nodes.remove(outgroup_name)
    kwargs['reduced_nodes']=reduced_nodes
    kwargs['sample_per_pop']=sample_per_pop
    kwargs['nreps']=nreps
    kwargs['theta']=theta
    kwargs['sites']=sites
    kwargs['recomb_rate']=recomb_rate 
    kwargs['ms_file']=ms_file
    kwargs['treemix_file']=treemix_file
    kwargs['blocksize_empirical_covariance']=blocksize_empirical_covariance
    kwargs['ms_variance_correction']=ms_variance_correction
    
    
    #makes a necessary transformation of the input(if the input is a filename or something).
    statistic=read_input(stages_to_go_through[0], input, full_nodes, reduced_nodes)
    
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
    
    return statistic

def write_output(stage, output):
    pass

def read_input(stage, input, full_nodes, reduced_nodes):
    if stage==1:
        return input
    if stage==2:
        return input
    if stage==3:
        return read_tree(input, reduced_nodes)
    if stage==4:
        return read_tree(input, full_nodes)
    if stage==5:
        return input
    if stage==6:
        return read_covariance_matrix(input, full_nodes)
    if stage==7:
        return read_covariance_matrix(input, reduced_nodes)
    if stage==8: #assuming that it comes in the pair (matrix or matrix filename, multiplier)
        return (read_covariance_matrix(input[0], nodes), input[1])
    assert False, 'The beginning state '+str(stage)+' is unknown.'
        
def read_covariance_matrix(input, nodes):
    if isinstance(input, basestring):
        return file_to_emp_cov(input, nodes)
    else:
        return input
    
def read_tree(input, nodes):
    if isinstance(input, basestring):
        if not ';' in input:
            input=read_one_line(filename=input)
            return identifier_to_tree_clean(input, nodes=nodes)
        else:
            return input
        
def read_one_line(filename):
    with open(filename, 'r') as f:
        return f.readline().rstrip()
    
if __name__=='__main__':
    print get_covariance(stages_to_go_through=[1,2,3,6], input=4, full_nodes=['s1','s2','s3','s4'])