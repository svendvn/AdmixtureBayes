from prior import prior
from likelihood import likelihood, n_mark
from scipy.stats import norm, multivariate_normal
from math import log



def initialize_posterior(emp_cov, M=None, use_skewed_distr=False, p=0.5, rescale=False):
    if M is None:
        M=n_mark(emp_cov)
    if rescale:
        emp_cov, multiplier = rescale_empirical_covariance(emp_cov)
        print 'multiplier is', multiplier
    def posterior(x,pks={}):
        #print tot_branch_length
        prior_value=prior(x,p=p, use_skewed_distr=use_skewed_distr,pks=pks)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        likelihood_value=likelihood(x, emp_cov,M=M)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value
    if rescale:
        return posterior, multiplier
    return posterior

def initialize_big_posterior(emp_cov, M=None, use_skewed_distr=False, p=0.5):
    if M is None:
        M=n_mark(emp_cov)
    def posterior(x,pks={}):
        #print tot_branch_length
        prior_value=prior(x,p=p, use_skewed_distr=use_skewed_distr,pks=pks)
        if prior_value==-float('inf'):
            return -float('inf'), prior_value
        likelihood_value=likelihood(x, emp_cov,M=M, pks=pks)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        prior_values=(pks['branch_prior'], pks['no_admix_prior'], pks['admix_prop_prior'], pks['top_prior'])
        covariance=pks['covariance']
        #pks['posterior']=prior_value+likelihood_value
        return likelihood_value, prior_value, prior_values, covariance
    return posterior
        

def initialize_prior_as_posterior():
    def posterior(x,pks={}):
        #print tot_branch_length
        prior_value=prior(x,pks=pks)
        if prior_value==-float('inf'):
            return prior_value
        pks['prior']=prior_value
        pks['likelihood']=0
        return 0,prior_value
    return posterior

def initialize_trivial_posterior():
    def posterior(x, pks={}):
        if isinstance(x, float):
            res= norm.logpdf(x)
            pks['prior']= res
            return res
#         elif isinstance(x, list) and isinstance(x[0], float):
#             res= multivariate_normal.logpdf(x)
#             pks['prior']= res
#             return res
        else:
            assert False, 'input in posterior was not recognizable.'
    return posterior

def print_pks(pks):
    for key, value in pks.items():
        print '\t',key,':',value

def call_post(posterior_function,tree, posterior_function_name='posterior', tree_name='tree'):
    pks={}
    print posterior_function_name+'('+tree_name+')=', posterior_function(tree, pks=pks)
    print_pks(pks)
    
def rescale_empirical_covariance(m):
    '''
    It is allowed to rescale the empirical covariance matrix such that the inferred covariance matrix takes values that are closer to the mean of the prior.
    '''
    
    n=m.shape[0]
    actual_trace=m.trace()
    expected_trace=log(n)/log(2)*n
    multiplier= expected_trace/actual_trace
    return m*multiplier, multiplier


if __name__=='__main__':
    import Rcatalogue_of_trees
    import Rtree_operations
    import Rtree_to_covariance_matrix
    
    true_tree=Rcatalogue_of_trees.tree_good
    ref_tree=Rtree_operations.create_trivial_tree(4, 0.2)
    nodes=Rtree_operations.get_trivial_nodes(4)
    
    true_cov=Rtree_to_covariance_matrix.make_covariance(true_tree, nodes)
    ref_cov=Rtree_to_covariance_matrix.make_covariance(ref_tree,nodes)
    
    true_posterior=initialize_posterior(true_cov)
    ref_posterior=initialize_posterior(ref_cov)

        
    
    call_post(true_posterior, ref_tree, 'true_posterior', 'ref_tree')
    call_post(true_posterior, true_tree, 'true_posterior', 'true_tree')
    call_post(ref_posterior, ref_tree, 'ref_posterior', 'ref_tree')
    call_post(ref_posterior, true_tree, 'ref_posterior', 'true_tree')
    
    print '''----fixed n'----- '''
    
    true_posterior=initialize_posterior(true_cov,M=10)
    ref_posterior=initialize_posterior(ref_cov,M=10)
    
    call_post(true_posterior, ref_tree, 'true_posterior', 'ref_tree')
    call_post(true_posterior, true_tree, 'true_posterior', 'true_tree')
    call_post(ref_posterior, ref_tree, 'ref_posterior', 'ref_tree')
    call_post(ref_posterior, true_tree, 'ref_posterior', 'true_tree')