from tree_to_covariance_matrix import make_covariance

from scipy.stats import wishart


def likelihood(tree, emp_cov):
    par_cov=make_covariance(tree)
    try:
        d=wishart.pdf(emp_cov, df=1, scale= par_cov)
    except ValueError:
        return 0
    return d
        
    

if __name__=="__main__":
    from tree_operations import 