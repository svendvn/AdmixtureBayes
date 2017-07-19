from Rtree_to_covariance_matrix import make_covariance

from scipy.stats import wishart
from numpy.linalg import eig, LinAlgError

def n_mark(cov_mat):
    m=cov_mat.shape[0]
    w,_=eig(cov_mat)
    S=sum(w)
    print S
    USS=sum(w*w)
    print USS
    return (m-1)*S**2 / ( (m+1)*USS-S**2 )
    

def likelihood(x, emp_cov, nodes=None, M=12, pks={}):
    tree, add= x
    r=emp_cov.shape[0]
    if nodes is None:
        nodes=["s"+str(i) for i in range(1,r+1)]
    par_cov=make_covariance(tree, nodes)
    #print (par_cov, emp_cov, add)
    pks['covariance']=par_cov
    if par_cov is None:
        print 'illegal tree'
        return -float('inf')
    try:
        #print emp_cov-add
        #print add
        #print par_cov
        d=wishart.logpdf(r*M*(emp_cov-add), df=r*M-1, scale= par_cov)
    except (ValueError, LinAlgError) as e:
        #print "illegal par_cov matrix or to large add"
        #print e
        return -float("inf")
    return d
        
    

if __name__=="__main__":
    import Rcatalogue_of_trees
    import Rtree_operations
    true_tree=Rcatalogue_of_trees.tree_good
    s_tree=Rtree_operations.create_trivial_tree(4)
    emp_cov=make_covariance(true_tree, Rtree_operations.get_trivial_nodes(4))
    nmark=n_mark(emp_cov)
    print 'm,nmark',emp_cov.shape[0], nmark
    print emp_cov
    print make_covariance(s_tree, Rtree_operations.get_trivial_nodes(4))
    print 'likelihood(M=12)', likelihood(s_tree, emp_cov)
    print 'likelihood(M=nmark)', likelihood(s_tree, emp_cov,M=nmark)
#     from tree_operations import make_flat_list_no_admix
#     from numpy import diag
#     N=5
#     nodes=["s"+str(i) for i in range(1,N+1)]
#     tree=make_flat_list_no_admix(5)
#     print tree
#     emp_cov=diag([0.5]*N)
#     print emp_cov
#     print likelihood(tree,emp_cov)