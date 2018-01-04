from generate_prior_covariance import generate_covariance
from covariance_estimator import initor
from covariance_EM import EmEstimator
from covariance_scaled import ScaledEstimator
from covariance_indirect_correction import IndirectEstimator
from covariance_simulation import Simulator

def create_initial_Sigma_generator(streng):
    if streng=='default':
        return fixed_initial_Sigma(None)
    elif streng=='random':
        return random_initial_Sigma(n, scale)

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
                   Simulator_from_file=''):

    if isinstance(initial_Sigma_generator, basestring):
        initial_Sigma_generator=create_initial_Sigma_generator(initial_Sigma_generator)
        

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
                 bias_c_weight='default')
    if indirect_correction:
        simulator=Simulator(ns,estimator=est, fixed_seed=Simulator_fixed_seed,  load_from_file=Simulator_from_file)
        est2=IndirectEstimator(no_its=Indirect_its,
                               s=Indirect_multiplier_s,
                               initial_Sigma_generator=initial_Sigma_generator,
                               estimator=est,
                               simulator=simulator)
    else:
        est2=est
        
    if no_repeats_of_cov_est>1:
        est3=RepeatEstimator(est2, no_repeats)
    else:
        est3=est2
    
    return est3