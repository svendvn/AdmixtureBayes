import subprocess
from numpy import array, mean, zeros, diag, sum, arcsin, sqrt
from reduce_covariance import reduce_covariance
#from optimize_empirical_matrix import full_maximization, transform_allele_freqs
from copy import deepcopy
import warnings
from covariance_EM import EmEstimator
from covariance_estimator import initor
from covariance_scaled import ScaledEstimator
from covariance_indirect_correction import IndirectEstimator











# def alleles_to_cov_optimizing(p, 
#                               names, 
#                               pop_sizes=None,
#                               nodes=None,
#                               cutoff=0.01):
#     p=reorder_rows(p, names, nodes)
#     print 'p',p
#     xs=p
#     ns=deepcopy(p)
#     
#     for n in range(len(pop_sizes)):
#         popsize=pop_sizes[n]
#         xs[n,:]*=popsize
#         print 'xs',xs
#         ns[n,:]=popsize
#         
#     print 'xs',xs
#     print 'ns',ns
#     
#     trans=transform_allele_freqs(cutoff)
#     
#     return full_maximization(xs,ns, trans=trans)
#     
#         
# 
# def alleles_to_cov(p,
#                    names, 
#                    pop_sizes=None, 
#                    reduce_method=['no', 'average', 'outgroup'], 
#                    variance_correction=['None','unbiased','mle'], 
#                    nodes=None, 
#                    arcsin_transform=False, 
#                    method_of_weighing_alleles=['None', 'Jade','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product','Jade-o'], 
#                    reducer='',
#                    jade_cutoff=1e-5,
#                    reduce_also=False,
#                    bias_c_weight='default'):
#     p=reorder_rows(p, names, nodes)
#     
#     if not isinstance(variance_correction, basestring):
#         type_of_scaling=variance_correction[0]
#     else:
#         type_of_scaling=variance_correction
#     
#     if not isinstance(reduce_method, basestring):
#         reduce_method=reduce_method[0]
#     if not isinstance(method_of_weighing_alleles, basestring):
#         method_of_weighing_alleles=method_of_weighing_alleles[0]
#     
#     if arcsin_transform:
#         p2=arcsin(sqrt(p))*2
#     else:
#         p2=p
#         
#     if reduce_method=='outgroup':
#         n_outgroup=next((n for n, e in enumerate(nodes) if e==reducer))
#         #print n_outgroup
#         p2=p2-p2[n_outgroup,:]
#     elif reduce_method=='average':
#         n_outgroup=next((n for n, e in enumerate(nodes) if e==reducer))
#         total_mean2=mean(p, axis=0)
#         if arcsin_transform:
#             total_mean2=arcsin(sqrt(total_mean2))*2
#         p2=p2-total_mean2
#     else:
#         n_outgroup=None
#     
#     if method_of_weighing_alleles=='Jade':
#         mu=mean(p, axis=0)
#         
#         i=array([v > jade_cutoff and v<1.0-jade_cutoff for v in mu ])
#         p2=p2[:,i]/sqrt(mu[i]*(1.0-mu[i]))
#         #p=p[:,i]
#     if method_of_weighing_alleles=='Jade-o':
#         mu=p[n_outgroup,:]
#         
#         i=array([v > jade_cutoff and v<1.0-jade_cutoff for v in mu ])
#         #p=p[:,i]
#         p2=p2[:,i]/sqrt(mu[i]*(1.0-mu[i]))
#         
#     m=p2.dot(p2.T)/p2.shape[1]
#     
#     if bias_c_weight=='default':
#         bias_c_weight=default_scale_dic[method_of_weighing_alleles]
#     
#     if method_of_weighing_alleles != 'None' and method_of_weighing_alleles!='Jade' and method_of_weighing_alleles != 'Jade-o':
#         m=m/m_scaler(m, method_of_weighing_alleles, p, n_outgroup)
#     
#     if reduce_also:
#         m=reduce_covariance(m, n_outgroup)
#         if type_of_scaling!='None':
#             m=other_bias_correction(m,p, pop_sizes,n_outgroup)
#     elif type_of_scaling!='None':
#         changer=bias_correction(m,p, pop_sizes,n_outgroup, type_of_scaling=type_of_scaling)/m_scaler(m, bias_c_weight, p, n_outgroup)
#         m=m/m_scaler(m, method_of_weighing_alleles, p, n_outgroup)
#         print 'm',reduce_covariance(m,n_outgroup)
#         print 'changer', reduce_covariance(changer, n_outgroup)
#         m-=changer
#         
#     return m



#wrapper to handle both treemix file input and numpy array input
def to_cov(snps='',
           names=None,
           reduce_method=['no', 'average', 'outgroup'], 
           variance_correction=['None', 'unbiased', 'mle'], 
           indirect_correction=False,
           nodes=None, 
           arcsin_transform=False, 
           method_of_weighing_alleles=['None', 'Jade','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product','Jade-o', 'EM'], 
           reducer='',
           jade_cutoff=1e-5,
           reduce_also=False,
           bias_c_weight='default'):
    
    if isinstance(snps, basestring):
        new_filename=make_uncompressed_copy(snps)
        xij, nij, names, pop_sizes, minors, total_sum= read_freqs(new_filename)
        allele_counts=array(allele_counts).T
    else:
        alle_counts=snps
    

        
        
    
    
    

def treemix_to_cov(filename='treemix_in.txt.gz', 
                   reduce_method=['no', 'average', 'outgroup'], 
                   reducer='', 
                   variance_correction=False, 
                   nodes=None,
                   arcsin_transform=False,
                   method_of_weighing_alleles='None',
                   jade_cutoff=1e-5
                   ):
#unzip
    new_filename=make_uncompressed_copy(filename)
#     print 'FILENAME', filename
#     print 'NEW FILENAME', new_filename
    
    allele_counts, names, pop_sizes, minors, total_sum= read_freqs(new_filename)
    
    p=array(allele_counts)
    p=p.T
    
    return alleles_to_cov(p, 
                          names, 
                          pop_sizes, 
                          reduce_method, 
                          variance_correction, 
                          nodes, 
                          arcsin_transform, 
                          method_of_weighing_alleles, 
                          reducer,
                          jade_cutoff)



if __name__=='__main__':
    
    if False:
        m=adjuster([0.1,0.2])
        assert sum(m-array([[0.075,-0.075],[-0.075,0.075]]), axis=None)==0, 'adjusting wrong'
        from numpy import set_printoptions
        set_printoptions(precision=5)
        filename='sletmig/_treemix_in.txt.gz'
        nodes=['s'+str(i) for i in range(1,10)]+['out']
        print treemix_to_cov(filename, reduce_method='outgroup', reducer='out', variance_correction=False, nodes=nodes, arcsin_transform=False, method_of_weighing_alleles='outgroup_product')
        
        from load_data import read_data
        
        print read_data(filename, blocksize=1, nodes=nodes, variance_correction=True, normalize=False, reduce_also=True, reducer='out', return_muhat=False)
    
    
    
    if True:
        import numpy as np

        from scipy.stats import uniform
        def is_pos_def(x):
            print np.linalg.eigvals(x)
            return np.all(np.linalg.eigvals(x) > 0)
        p=array([uniform.rvs(size=1000) for _ in xrange(5)])
        a=-other_bias_correction(np.zeros((4,4)), p, [5,5,5,5,5],0)
        print a
        print is_pos_def(a)