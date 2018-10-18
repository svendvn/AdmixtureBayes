from generate_prior_covariance import generate_covariance
from covariance_estimator import initor
from covariance_EM import EmEstimator
from covariance_scaled import ScaledEstimator
from covariance_indirect_correction import IndirectEstimator
from covariance_simulation import Simulator
from covariance_estimator import RepeatEstimator

from Rtree_to_covariance_matrix import make_covariance


def create_initial_Sigma_generator(n, streng):
    key=streng.keys()[0]
    if key=='default':
        return fixed_initial_Sigma(None)
    elif key=='random':
        return random_initial_Sigma(n)
    elif key=='start':
        print streng[key]
        cov=make_covariance(streng[key][0][0], node_keys=streng[key][1])+streng[key][0][1]
        return fixed_initial_Sigma(cov)
    

def fixed_initial_Sigma(initial_Sigma):
    def f():
        return initial_Sigma
    return f

def random_initial_Sigma(n, scale='beta'):
    def f():
        return generate_covariance(n, scale)
    return f
    
# def make_estimator(reduce_method=['no', 'average', 'outgroup'], 
#                    variance_correction=['None', 'unbiased', 'mle'], 
#                    indirect_correction=False,
#                    nodes=None, 
#                    arcsin_transform=False, 
#                    method_of_weighing_alleles=['None', 'Jade','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product','Jade-o', 'EM'], 
#                    reducer='',
#                    jade_cutoff=1e-5,
#                    reduce_also=True,
#                    bias_c_weight='default',
#                    EM_maxits=100,
#                    Indirect_its=100,
#                    EM_alpha=1.0,
#                    Indirect_multiplier_s=1,
#                    Simulator_fixed_seed=True,
#                    Simulator_from_file='',
#                    initial_Sigma_generator=None,
#                    no_repeats_of_cov_est=1,
#                    ns=None, #only necessary if indirect estimation is used.
#                    ):
    
def make_estimator(reduce_method, 
                   variance_correction, 
                   indirect_correction,
                   nodes, 
                   arcsin_transform, 
                   method_of_weighing_alleles, 
                   reducer,
                   jade_cutoff,
                   reduce_also,
                   bias_c_weight,
                   EM_maxits,
                   Indirect_its,
                   EM_alpha,
                   Indirect_multiplier_s,
                   Simulator_fixed_seed,
                   initial_Sigma_generator,
                   no_repeats_of_cov_est,
                   ns, #only necessary if indirect estimation is used.
                   Simulator_from_file='',
                   locus_filter_on_simulated=None,
                   add_variance_correction_to_graph=False,
                   prefix='',
                   save_variance_correction=True):

    n=len(nodes)-int(reduce_also)
    initial_Sigma_generator=create_initial_Sigma_generator(n, initial_Sigma_generator)
        

    variance_correction=initor(variance_correction)
    
    method_of_weighing_alleles=initor(method_of_weighing_alleles)
    
    if method_of_weighing_alleles=='EM':
        est=EmEstimator(maxiter=EM_maxits,
                        alpha=EM_alpha,
                        initial_Sigma_generator=initial_Sigma_generator)
    else:
        est=ScaledEstimator(reduce_method=reduce_method,
                 scaling=method_of_weighing_alleles,
                 reduce_also=reduce_also,
                 variance_correction=variance_correction,
                 jade_cutoff=1e-5,
                 bias_c_weight='default',
                 add_variance_correction_to_graph=add_variance_correction_to_graph,
                 prefix_for_saving_variance_correction=prefix,
                 save_variance_correction=save_variance_correction)
    if indirect_correction:
        simulator=Simulator(ns, multiplier=Indirect_multiplier_s, estimator=est, fixed_seed=Simulator_fixed_seed,  load_from_file=Simulator_from_file, locus_filter=locus_filter_on_simulated)
        print 'MULTIPLIER', Indirect_multiplier_s
        est2=IndirectEstimator(no_its=Indirect_its,
                               s=Indirect_multiplier_s,
                               initial_Sigma_generator=initial_Sigma_generator,
                               estimator=est,
                               simulator=simulator)
    else:
        est2=est
        
    if no_repeats_of_cov_est>1:
        est3=RepeatEstimator(est2, no_repeats_of_cov_est)
    else:
        est3=est2
    
    return est3

if __name__=='__main__':
    from scipy.stats import binom
    import numpy as np
    pop_sizes=np.ones((5,4))*3
    xs=binom.rvs(3,0.5, size=pop_sizes.shape)
    from covariance_indirect_correction import search_sigmas, search_choleskys
    from covariance_simulation import Simulator
    from tree_to_data import get_xs_and_ns_from_treemix_file
    from construct_filter_choices import make_filter
    from covariance_simulation import Simulator
    import numpy as np
    treemix_file='../../../../Dropbox/Bioinformatik/AdmixtureBayes/likelihood_tests/tmp_files/_treemix_in.txt.gz'
    locus_filter=make_filter(outgroup_name='out', filter_type='all_pops')
    outp=get_xs_and_ns_from_treemix_file(treemix_file, locus_filter)
    locus_filter2=make_filter(filter_type='all_pops')
    xs,ns,_=outp
    print xs.shape
    print ns.shape
    est=make_estimator(reduce_method='outgroup', 
                       variance_correction='None', 
                       indirect_correction=False,
                       nodes=['s1','s2','s3','out'], 
                       arcsin_transform=False, 
                       method_of_weighing_alleles='average_sum', 
                       reducer='out',
                       jade_cutoff=1e-5,
                       reduce_also=True,
                       bias_c_weight='average_sum',
                       EM_maxits=1,
                       Indirect_its=1,
                       EM_alpha=1,
                       Indirect_multiplier_s=1,
                       Simulator_fixed_seed=False,
                       initial_Sigma_generator={'default':''},
                       no_repeats_of_cov_est=1,
                       ns=np.array(pop_sizes), #only necessary if indirect estimation is used.
                       Simulator_from_file='',
                       locus_filter_on_simulated=None)
        
    print est(xs,ns)
    
    