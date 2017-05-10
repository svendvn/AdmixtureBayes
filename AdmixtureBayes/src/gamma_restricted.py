from scipy.stats import gamma

def conditional_gamma_rvs(mean,variance, upper_limit):

def rvs(t,delta_L, variance=0.1):
    if delta_L<0:
        mean=-delta_L
        return t-conditional_gamma(mean, variance=0.1, upper_limit)
    else:
        return t+gamma.rvs()
    
def clean_gamma_rvs(mean, variance):
    
        
    