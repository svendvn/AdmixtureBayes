from prior import prior
from likelihood import likelihood
from scipy.stats import norm, multivariate_normal

def initialize_posterior(emp_cov):
    def posterior(x,pks={}):
        #print tot_branch_length
        prior_value=prior(x,pks=pks)
        if prior_value==-float('inf'):
            return prior_value
        likelihood_value=likelihood(x, emp_cov)
        pks['prior']=prior_value
        pks['likelihood']=likelihood_value
        return prior_value+likelihood_value
    return posterior

def initialize_prior_as_posterior():
    def posterior(x,pks={}):
        #print tot_branch_length
        prior_value=prior(x,pks=pks)
        if prior_value==-float('inf'):
            return prior_value
        pks['prior']=prior_value
        pks['likelihood']=''
        return prior_value
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