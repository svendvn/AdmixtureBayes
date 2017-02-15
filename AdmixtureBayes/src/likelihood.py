from tree_to_covariance_matrix import make_covariance

from scipy.stats import wishart


def likelihood(tree, emp_cov, nodes=None):
    r=emp_cov.shape[0]
    if nodes is None:
        nodes=["s"+str(i) for i in range(1,r+1)]
    print nodes
    par_cov=make_covariance(tree, nodes)
    print par_cov
    print tree
    try:
        d=wishart.pdf(emp_cov, df=r, scale= par_cov)
    except ValueError:
        print "illegal par_cov"
        return 0
    return d
        
    

if __name__=="__main__":
    from tree_operations import make_flat_list_no_admix
    from numpy import diag
    N=5
    nodes=["s"+str(i) for i in range(1,N+1)]
    tree=make_flat_list_no_admix(5)
    print tree
    emp_cov=diag([0.5]*N)
    print emp_cov
    print likelihood(tree,emp_cov)