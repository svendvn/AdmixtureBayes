from tree_statistics import identifier_to_tree_clean, unique_identifier_and_branch_lengths
from tree_to_data import file_to_emp_cov, reduce_covariance, ms_to_treemix3, call_ms_string, tree_to_ms_command, emp_cov_to_file
from generate_prior_trees import simulate_number_of_admixture_events, generate_phylogeny
from Rtree_operations import add_outgroup, get_number_of_leaves, scale_tree
from scipy.stats import expon
from Rtree_to_covariance_matrix import make_covariance
from load_data import read_data
from copy import deepcopy
from math import log





def simulate_tree_wrapper(nk_tuple, **kwargs):
    return generate_phylogeny(size= nk_tuple[0], 
                              admixes=nk_tuple[1], 
                              leaf_nodes= kwargs['before_added_outgroup_nodes'], 
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
        covariance=add_wishart_noise(covariance, kwargs['df_of_wishart_noise_to_covariance'])
    return covariance
        
def add_wishart_noise(matrix, df):
    r=matrix.shape[0]
    m=wishart.rvs(df=r*df-1, scale=matrix/(r*df))
    return m

def reduce_covariance_wrapper(full_covariance, **kwargs):
    reduce_node_index=next((i for i,s in enumerate(kwargs['full_nodes']) if s==kwargs['outgroup_name']))
    return reduce_covariance(full_covariance, reduce_node_index)
    
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
                   filename2=kwargs['treemix_file'],
                   nodes=kwargs['full_nodes'])
    return filename_gz

def empirical_covariance_wrapper(snp_data_file, **kwargs):
    cov=read_data(snp_data_file, 
                     outgroup= '',
                     blocksize=kwargs['blocksize_empirical_covariance'],
                     nodes=kwargs['full_nodes'], 
                     noss=kwargs['ms_variance_correction'])
    return cov
    
def scale_tree_wrapper(tree, **kwargs):
    return scale_tree(tree, kwargs['scale_tree_factor'])

def normaliser_wrapper(covariance, **kwargs):
    return rescale_empirical_covariance(covariance)




dictionary_of_transformations={
    (1,2):simulate_number_of_admixture_events_wrapper,
    (2,3):simulate_tree_wrapper,
    (3,4):add_outgroup_wrapper,
    (4,5):scale_tree_wrapper,
    (3,5):scale_tree_wrapper,
    (3,7):theoretical_covariance_wrapper,
    (4,7):theoretical_covariance_wrapper,
    (5,7):theoretical_covariance_wrapper,
    (3,6):ms_simulate_wrapper,
    (4,6):ms_simulate_wrapper,
    (5,6):ms_simulate_wrapper,
    (6,7):empirical_covariance_wrapper,
    (7,8):reduce_covariance_wrapper,
    (7,9):normaliser_wrapper,
    (8,9):normaliser_wrapper
    }

dictionary_of_reasonable_names={
    1:'number_of_leaves',
    2:'leaves_admixtures',
    3:'true_tree',
    4:'true_tree_with_outgroup',
    5:'scaled_true_tree',
    6:'SNP_data',
    7:'covariance',
    8:'covariance_without_reduce_name',
    9:'covariance_and_multiplier'}

def write_one_line_to_file(filename, value):
    with open(filename,'w') as f:
        f.write(value)

def save_stage(value, stage_number, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes):
    save_word=dictionary_of_reasonable_names[stage_number]
    filename=prefix+save_word+'.txt'
    if stage_number==1:
        write_one_line_to_file(filename, str(value))
    elif stage_number==2:
        write_one_line_to_file(filename, str(value))
    elif stage_number==3:
        write_one_line_to_file(filename, unique_identifier_and_branch_lengths(value))
    elif stage_number==4:
        write_one_line_to_file(filename, unique_identifier_and_branch_lengths(value))
    elif stage_number==5:
        write_one_line_to_file(filename, unique_identifier_and_branch_lengths(value))
    elif stage_number==6:
        print 'file is already made elsewhere'
    elif stage_number==7:
        emp_cov_to_file(value, filename, full_nodes)
    elif stage_number==8:
        emp_cov_to_file(value, filename, after_reduce_nodes)
    else:
        emp_cov_to_file(value[0], filename, after_reduce_nodes)
        with open(filename, 'a') as f:
            f.write('multiplier='+str(value[1]))

def rescale_empirical_covariance(m):
    '''
    It is allowed to rescale the empirical covariance matrix such that the inferred covariance matrix takes values that are closer to the mean of the prior.
    '''
    
    n=m.shape[0]
    actual_trace=m.trace()
    expected_trace=log(n)/log(2)*n
    multiplier= expected_trace/actual_trace
    return m*multiplier, multiplier


def get_covariance(stages_to_go_through, input, full_nodes=None,
                   skewed_admixture_prior_sim=False,
                   p=0.5,
                   outgroup_name=None,
                   add_wishart_noise_to_covariance=False,
                   df_of_wishart_noise_to_covariance=1000,
                   reduce_covariance_node=None,
                   sample_per_pop=50, nreps=2, 
                   theta=0.4, sites=500000, recomb_rate=1,
                   ms_file=None,
                   treemix_file=None,
                   blocksize_empirical_covariance=100,
                   ms_variance_correction=False,
                   scale_tree_factor=0.05,
                   save_stages=range(1,6)+range(7,10),
                   prefix='tmp'):
    
    if prefix[-1]!='_':
        prefix+='_'
    
    if ms_file is None:
        ms_file=prefix+'ms.txt'
    
    if treemix_file is None:
        treemix_file=prefix+'treemix'
    
    kwargs={}
    kwargs['skewed_admixture_prior_sim']=skewed_admixture_prior_sim
    kwargs['p']=p
    kwargs['outgroup_name']=outgroup_name
    kwargs['add_wishart_noise_to_covariance']=add_wishart_noise_to_covariance
    kwargs['df_of_wishart_noise_to_covariance']=df_of_wishart_noise_to_covariance
    kwargs['full_nodes']=full_nodes
    before_added_outgroup_nodes=deepcopy(full_nodes)
    after_reduce_nodes=deepcopy(full_nodes)
    if outgroup_name is not None and outgroup_name in before_added_outgroup_nodes:
        before_added_outgroup_nodes.remove(outgroup_name)
    if reduce_covariance_node is not None and reduce_covariance_node in after_reduce_nodes:
        after_reduce_nodes.remove(reduce_covariance_node)
    kwargs['after_reduce_nodes']=after_reduce_nodes
    kwargs['before_added_outgroup_nodes']=before_added_outgroup_nodes
    kwargs['sample_per_pop']=sample_per_pop
    kwargs['nreps']=nreps
    kwargs['theta']=theta
    kwargs['sites']=sites
    kwargs['recomb_rate']=recomb_rate 
    kwargs['ms_file']=ms_file
    kwargs['treemix_file']=treemix_file
    kwargs['blocksize_empirical_covariance']=blocksize_empirical_covariance
    kwargs['ms_variance_correction']=ms_variance_correction
    kwargs['scale_tree_factor']=scale_tree_factor
    
    
    #makes a necessary transformation of the input(if the input is a filename or something).
    statistic=read_input(stages_to_go_through[0], input, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
    
    if stages_to_go_through[0] in save_stages:
        save_stage(statistic, stages_to_go_through[0], prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
    
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        print (stage_from, stage_to)
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
        if stage_to in save_stages:
            save_stage(statistic, stage_to, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
        print statistic
    
    return statistic

def write_output(stage, output):
    pass

def read_input(stage, input, full_nodes, before_added_outgroup_nodes, after_reduce_nodes):
    if stage==1:
        return int(input)
    if stage==2:
        return map(int, input[1:-1].split(','))
    if stage==3:
        return read_tree(input, before_added_outgroup_nodes)
    if stage==4:
        return read_tree(input, full_nodes)
    if stage==5:
        return read_tree(input, full_nodes)
    if stage==6:
        return input
    if stage==7:
        return read_covariance_matrix(input, full_nodes)
    if stage==8:
        return read_covariance_matrix(input, after_reduce_nodes)
    if stage==9 : #assuming that it comes in the pair (matrix or matrix filename, multiplier)
        return (read_covariance_matrix(input[0], full_nodes), input[1])
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
            return identifier_to_tree_clean(input, leaves=nodes)
        else:
            return identifier_to_tree_clean(input, leaves=nodes)
    else:
        return input
        
def read_one_line(filename):
    with open(filename, 'r') as f:
        return f.readline().rstrip()
    
if __name__=='__main__':
    print get_covariance(stages_to_go_through=[1,2,3,4,5,6,7,8,9], input=4, full_nodes=['s1','s2','s3','s4','outgroup_name'])