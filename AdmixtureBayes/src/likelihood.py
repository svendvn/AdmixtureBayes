from Rtree_to_covariance_matrix import make_covariance

from scipy.stats import wishart


def likelihood(tree, emp_cov, nodes=None, M=12):
    r=emp_cov.shape[0]
    if nodes is None:
        nodes=["s"+str(i) for i in range(1,r+1)]
    par_cov=make_covariance(tree, nodes)
    if par_cov is None:
        return -float('inf')
    try:
        d=wishart.logpdf(r*M*emp_cov, df=r*M-1, scale= par_cov)
    except ValueError:
        print "illegal par_cov matrix"
        return -float("inf")
    return d
        
    

if __name__=="__main__":
    import Rcatalogue_of_trees
    import Rtree_operations
    true_tree=Rcatalogue_of_trees.tree_good
    s_tree=Rtree_operations.create_trivial_tree(4)
    emp_cov=make_covariance(true_tree, Rtree_operations.get_trivial_nodes(4))
    print emp_cov
    print make_covariance(s_tree, Rtree_operations.get_trivial_nodes(4))
    print likelihood(s_tree, emp_cov)
#     from tree_operations import make_flat_list_no_admix
#     from numpy import diag
#     N=5
#     nodes=["s"+str(i) for i in range(1,N+1)]
#     tree=make_flat_list_no_admix(5)
#     print tree
#     emp_cov=diag([0.5]*N)
#     print emp_cov
#     print likelihood(tree,emp_cov)