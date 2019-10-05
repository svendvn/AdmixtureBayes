from tree_statistics import identifier_to_tree_clean, unique_identifier_and_branch_lengths, generate_predefined_list_string
from tree_to_data import (file_to_emp_cov, reduce_covariance, ms_to_treemix3, call_ms_string, 
                          tree_to_ms_command, emp_cov_to_file, time_adjusted_tree_to_ms_command,
                           get_xs_and_ns_from_freqs, 
                          get_xs_and_ns_from_treemix_file, order_covariance, reorder_covariance, reorder_reduced_covariance)
from generate_prior_trees import simulate_number_of_admixture_events, generate_phylogeny
from generate_sadmix_trees import generate_sadmix_tree
from Rtree_operations import add_outgroup, get_number_of_leaves, scale_tree, time_adjust_tree
from scipy.stats import expon, wishart
from Rtree_to_covariance_matrix import make_covariance
from load_data import read_data
from copy import deepcopy
from math import log
from rescale_covariance import rescale_empirical_covariance
import time
from brownian_motion_generation import produce_p_matrix, calculate_covariance_matrix_from_p, simulate_with_binomial, remove_non_snps
from numpy import array, ones, loadtxt, savetxt
from construct_estimator_choices import make_estimator
import warnings





def simulate_tree_wrapper(nk_tuple, **kwargs):
    if kwargs['sadmix']:
        return generate_sadmix_tree(nk_tuple[0], no_sadmixes=nk_tuple[1], nodes=kwargs['before_added_outgroup_nodes'], starting_admixes=0)
    return generate_phylogeny(size= nk_tuple[0], 
                              admixes=nk_tuple[1], 
                              leaf_nodes= kwargs['before_added_outgroup_nodes'], 
                              skewed_admixture_prior=kwargs['skewed_admixture_prior_sim'])
    
def add_outgroup_wrapper(tree_without_outgroup, **kwargs):
    v=expon.rvs()
    new_length_1, new_length_2= v/2, v/2
    with open(kwargs['add_file'], 'w') as f:
        f.write(str(v))
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

def empirical_covariance_wrapper(snp_data_file, **kwargs):
    xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file, kwargs['locus_filter'])
    return xnn_to_covariance_wrapper(xnn_tuple, **kwargs)

def alleles_to_cov_wrapper(ps, **kwargs):
    xnn_tuple=get_xs_and_ns_from_freqs(ps, kwargs['sample_per_pop'], kwargs['locus_filter'])
    return xnn_to_covariance_wrapper(xnn_tuple, **kwargs)

def empirical_covariance_wrapper_directly(snp_data_file, **kwargs):
    xnn_tuple=get_xs_and_ns_from_treemix_file(snp_data_file, kwargs['locus_filter'])
    return xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs)

def alleles_to_cov_wrapper_directly(ps, **kwargs):
    xnn_tuple=get_xs_and_ns_from_freqs(ps, kwargs['sample_per_pop'], kwargs['locus_filter'])
    return xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs)

def xnn_to_covariance_wrapper(xnn_tuple, **kwargs):
    est_args=kwargs['est']
    xnn_tuple=order_covariance(xnn_tuple, outgroup=est_args['reducer'])
    xs,ns,names=xnn_tuple
   # for e,a in est_args.items():
     #   print e,a
    est= make_estimator(reduce_method='outgroup', 
                   reduce_also=False,
                   ns=ns,**est_args)
    cov=est(xs,ns)
    cov=reorder_covariance(cov, names, kwargs['full_nodes'])
    if ('add_variance_correction_to_graph' in est_args and 
        est_args['add_variance_correction_to_graph']):
        filename=est_args['prefix']+'variance_correction.txt'
        #print 'CHaningin VC'
        vc=loadtxt(filename)
        vc=reorder_covariance(vc, names, kwargs['full_nodes'])
        savetxt(filename, vc)
        
        
    return cov
    

def xnn_to_covariance_wrapper_directly(xnn_tuple, **kwargs):
    est_args=kwargs['est']
    #print 'kwargs:'
    #for i,j in est_args.items():
    #    print i,j
    #print xnn_tuple
    #print xnn_tuple
    #print 'shapes',(xnn_tuple[0].shape,xnn_tuple[1].shape)
    xnn_tuple=order_covariance(xnn_tuple, outgroup=est_args['reducer'])
    xs,ns,names=xnn_tuple
    

    est= make_estimator(reduce_method='outgroup', 
                   reduce_also=True,
                   ns=ns,**est_args)
    extra_info_dic={}
    cov=est(xs,ns, extra_info_dic)
    cov=reorder_reduced_covariance(cov, names, est_args['nodes'], outgroup=est_args['reducer'])
    if ('add_variance_correction_to_graph' in est_args and 
        est_args['add_variance_correction_to_graph'] and
        'save_variance_correction' in est_args and
        est_args['save_variance_correction']):
       # print 'CHaninging VC'
        filename=est_args['prefix']+'variance_correction.txt'
        vc=loadtxt(filename)
        vc=reorder_reduced_covariance(vc, names, est_args['nodes'], outgroup=est_args['reducer'])
        savetxt(filename, vc)
    if 'm_scale' in extra_info_dic:
        if 'mscale_file' in kwargs: 
           # print kwargs['mscale_file']
            with open(kwargs['mscale_file'], 'w') as f:
                txt=str(extra_info_dic['m_scale'])
                assert '\0' not in txt, 'binary content in m_scale file'
                f.write(txt)
        if 'return_also_mscale' in kwargs and kwargs['return_also_mscale']:
            return cov, extra_info_dic['m_scale']
    return cov
    

def add_wishart_noise(matrix, df):
    r=matrix.shape[0]
    m=wishart.rvs(df=df, scale=matrix/df)
    return m

def reduce_covariance_wrapper(full_covariance, **kwargs):
    reduce_node_index=next((i for i,s in enumerate(kwargs['full_nodes']) if s==kwargs['reduce_covariance_node']))
    return reduce_covariance(full_covariance, reduce_node_index)
    
def ms_simulate_wrapper(tree, **kwargs):
    no_pops= get_number_of_leaves(tree)
    #ms_command, minmaxv=tree_to_ms_command(tree,  #TO CHANGE BACK
    if kwargs['time_adjust']:
        ms_command=time_adjusted_tree_to_ms_command(tree,  #TO CHANGE BACK
                       sample_per_pop=kwargs['sample_per_pop'],
                       nreps=kwargs['nreps'],
                       theta=kwargs['theta'],
                       sites=kwargs['sites'],
                       recomb_rate=kwargs['recomb_rate'],
                       leaf_keys=kwargs['full_nodes'],
                       final_pop_size=kwargs['final_pop_size'])
    else:
        ms_command=tree_to_ms_command(tree,  #TO CHANGE BACK
                       sample_per_pop=kwargs['sample_per_pop'],
                       nreps=kwargs['nreps'],
                       theta=kwargs['theta'],
                       sites=kwargs['sites'],
                       recomb_rate=kwargs['recomb_rate'],
                       leaf_keys=kwargs['full_nodes'])
    #kwargs['pks']['minmaxv']=minmaxv  #TO CHANGE BACK
    #print ms_command
    call_ms_string(ms_command, kwargs['ms_file'])
    filename_gz=ms_to_treemix3(kwargs['ms_file'], 
                    samples_per_pop=kwargs['sample_per_pop'], 
                    no_pops=no_pops, 
                    n_reps=kwargs['nreps'], 
                    filename2=kwargs['treemix_file'],
                    nodes=kwargs['full_nodes'])
    return filename_gz



def simulate_brownian_motion_wrapper(tree, **kwargs):
    N=kwargs['nreps'] #number of independent SNPs. 
    ps=produce_p_matrix(tree,N, clip=(not kwargs['unbounded_brownian']), middle_start=kwargs['favorable_init_brownian'])
    return ps

def simulate_binomial_wrapper(p_dic, **kwargs):
    npop=kwargs['sample_per_pop']
    ps=simulate_with_binomial(p_dic, npop, p_clip_value=0)
    return ps

    

# def empirical_covariance_wrapper(snp_data_file, **kwargs):
#     if not kwargs['via_treemix']:
#         outgroup_number=any((n for n,e in enumerate(kwargs['full_nodes']) if e==kwargs['reduce_covariance_node']))
#         return calculate_covariance_matrix2(snp_data_file, 
#                                             samples_per_pop=kwargs['sample_per_pop'], 
#                                             no_pops=len(kwargs['full_nodes']), 
#                                             n_reps=kwargs['nreps'],
#                                             outgroup_number= outgroup_number)
#     cov=read_data(snp_data_file, 
#                      blocksize=kwargs['blocksize_empirical_covariance'],
#                      nodes=kwargs['full_nodes'], 
#                      noss=kwargs['ms_variance_correction'],
#                      outfile=kwargs['treemix_out_files'])
#     return cov

def remove_non_snps_wrapper(ps, **kwargs):
    outg=''
    if kwargs['filter_on_outgroup']:
        outg=kwargs['reduce_covariance_node']
    return remove_non_snps(ps, outg)
    
def scale_tree_wrapper(tree, **kwargs):
    if kwargs['time_adjust']:
        tree=time_adjust_tree(tree) 
    return scale_tree(tree, kwargs['scale_tree_factor'])

def normaliser_wrapper(covariance, **kwargs):
    return rescale_empirical_covariance(covariance, normalizer=kwargs['scale_goal'],
                                        outgroup_rate=kwargs['outgroup_branch_prior_rate'])




dictionary_of_transformations={
    (1,2):simulate_number_of_admixture_events_wrapper,
    (2,3):simulate_tree_wrapper,
    (3,4):add_outgroup_wrapper,
    (4,5):scale_tree_wrapper,
    (3,5):scale_tree_wrapper,
    (3,7):theoretical_covariance_wrapper,
    (4,7):theoretical_covariance_wrapper,
    (5,7):theoretical_covariance_wrapper,
    (3,21):simulate_brownian_motion_wrapper,
    (4,21):simulate_brownian_motion_wrapper,
    (5,21):simulate_brownian_motion_wrapper,
    (21,22):simulate_binomial_wrapper,
    (21,23):remove_non_snps_wrapper,
    (22,23):remove_non_snps_wrapper,
    (21,7):alleles_to_cov_wrapper,
    (21,8):alleles_to_cov_wrapper_directly,
    (22,7): alleles_to_cov_wrapper,
    (22,8): alleles_to_cov_wrapper_directly,
    (23,7): alleles_to_cov_wrapper,
    (23,8): alleles_to_cov_wrapper_directly,
    (3,6):ms_simulate_wrapper,
    (4,6):ms_simulate_wrapper,
    (5,6):ms_simulate_wrapper,
    (6,7):empirical_covariance_wrapper,
    (6,8):empirical_covariance_wrapper_directly,
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
        
def write_two_lines_to_file(filename, value1, value2):
    with open(filename, 'w') as f:
        f.write(value1+'\n'+value2)

def save_stage(value, stage_number, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes, filename=None):
    if filename is None:
        save_word=dictionary_of_reasonable_names[stage_number]
        filename=prefix+save_word+'.txt'
    if stage_number==1:
        write_one_line_to_file(filename, str(value))
    elif stage_number==2:
        write_one_line_to_file(filename, str(value))
    elif stage_number==3:
        write_two_lines_to_file(filename, ' '.join(before_added_outgroup_nodes), unique_identifier_and_branch_lengths(value, before_added_outgroup_nodes))
    elif stage_number==4:
        write_two_lines_to_file(filename, ' '.join(full_nodes), unique_identifier_and_branch_lengths(value, full_nodes))
    elif stage_number==5:
       # print full_nodes
        write_two_lines_to_file(filename, ' '.join(full_nodes), unique_identifier_and_branch_lengths(value, full_nodes))
    elif stage_number==6:
        pass
        #print 'file is already made elsewhere'
    elif stage_number==7:
        emp_cov_to_file(value, filename, full_nodes)
    elif stage_number==8:
        emp_cov_to_file(value, filename, after_reduce_nodes)
    elif stage_number in [21,22,23]:
        print 'Stage not saved'
    else:
        emp_cov_to_file(value[0], filename, after_reduce_nodes)
        with open(filename, 'a') as f:
            f.write('multiplier='+str(value[1]))




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
                   scale_tree_factor=0.05,
                   save_stages=range(1,6)+range(7,10),
                   prefix='tmp',
                   t_adjust_tree=False,
                   final_pop_size=100.0,
                   via_treemix=True,
                   sadmix=False,
                   scale_goal='min',
                   favorable_init_brownian=False,
                   unbounded_brownian=False,
                   filter_on_outgroup=False,
                   locus_filter=None,
                   outgroup_branch_prior_rate=1,
                   estimator_arguments={}, 
                   verbose_level='normal'):
    
    #if prefix[-1]!='_':
    #    prefix+='_'
    
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
    kwargs['sadmix']=sadmix
    before_added_outgroup_nodes=deepcopy(full_nodes)
    after_reduce_nodes=deepcopy(full_nodes)
    if outgroup_name is not None and outgroup_name in before_added_outgroup_nodes:
        before_added_outgroup_nodes.remove(outgroup_name)
    if reduce_covariance_node is not None and reduce_covariance_node in after_reduce_nodes:
        after_reduce_nodes.remove(reduce_covariance_node)
    kwargs['reduce_covariance_node']=reduce_covariance_node
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
    kwargs['scale_tree_factor']=scale_tree_factor
    kwargs['pks']={}
    kwargs['time_adjust']=t_adjust_tree
    kwargs['final_pop_size']=final_pop_size
    kwargs['via_treemix']=via_treemix
    kwargs['add_file']=prefix+'true_add.txt'
    kwargs['mscale_file']=prefix+'m_scale.txt'
    kwargs['scale_goal']=scale_goal
    kwargs['favorable_init_brownian']=favorable_init_brownian
    kwargs['unbounded_brownian']=unbounded_brownian
    kwargs['filter_on_outgroup']=filter_on_outgroup
    kwargs['outgroup_branch_prior_rate']=outgroup_branch_prior_rate
    kwargs['est']=estimator_arguments
    kwargs['locus_filter']=locus_filter
    
    start=time.time()
    #makes a necessary transformation of the input(if the input is a filename or something).
    statistic=read_input(stages_to_go_through[0], input, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
    
    if stages_to_go_through[0] in save_stages:
        save_stage(statistic, stages_to_go_through[0], prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
    
    for stage_from, stage_to in zip(stages_to_go_through[:-1], stages_to_go_through[1:]):
        #print (stage_from, stage_to)
        transformer_function=dictionary_of_transformations[(stage_from, stage_to)]
        statistic=transformer_function(statistic, **kwargs)
        if stage_to in save_stages:
            save_stage(statistic, stage_to, prefix, full_nodes, before_added_outgroup_nodes, after_reduce_nodes)
        #print statistic
    
    end=time.time()
    
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
        return (read_covariance_matrix(input, after_reduce_nodes), read_multiplier(input))
    assert False, 'The beginning state '+str(stage)+' is unknown.'
        
def read_multiplier(input):
    with open(input, 'r') as f:
        last_line=f.readlines()[-1]
        return float(last_line.split("=")[1])
        
def read_covariance_matrix(input, nodes):
    if isinstance(input, basestring):
        return file_to_emp_cov(input, nodes=nodes)
    else:
        return input
    
def read_tree(input, nodes):
    if isinstance(input, basestring):
        if not ';' in input:
            input=read_one_line_skip(filename=input)
            return identifier_to_tree_clean(input, leaves=generate_predefined_list_string(deepcopy(nodes)))
        else:
            return identifier_to_tree_clean(input, leaves=generate_predefined_list_string(deepcopy(nodes)))
    else:
        return input

def read_one_line_skip(filename):
    with open(filename, 'r') as f:
        lines=f.readlines()
        if len(lines[-1])>3:
            return lines[-1].rstrip()
        else:
            return lines[-2].rstrip()

def read_one_line(filename):
    with open(filename, 'r') as f:
        res=f.readline().rstrip()
    return res
    
if __name__=='__main__':
    #
    for _ in xrange(1):
        cov=get_covariance(stages_to_go_through=[2,3,4,5,6,7,8,9], input='(6,2)', full_nodes=['s1','s2','s3','s4','s5','s6','outgroup_name'],
                             outgroup_name='outgroup_name', reduce_covariance_node='outgroup_name',
                             nreps=500, t_adjust_tree=True, scale_tree_factor=0.01, ms_variance_correction=True, via_treemix=True,
                             add_wishart_noise_to_covariance=False)
        #print 'levels',levels
        import post_analysis
        other= post_analysis.get_true_posterior(outgroup='outgroup_name')
        resl=list(other)
        with open('resfile.txt', 'a') as f:
            f.write(','.join(map(str,resl))+'\n')
            
    
    
    
    
