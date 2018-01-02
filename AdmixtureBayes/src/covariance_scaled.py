from covariance_estimator import Estimator, initor

default_scale_dic={'None':'None',
                   'Jade-o':'outgroup_sum', 
                   'Jade':'average_sum',
                   'outgroup_sum':'outgroup_sum',
                   'outgroup_product':'outgroup_product',
                   'average_sum':'average_sum',
                   'average_product':'average_product'}

def m_scaler(m, scale_type, allele_freqs, n_outgroup=None):
    if scale_type=='None' or scale_type=='Jade' or scale_type=='Jade-o':
        return 1.0
    if scale_type.startswith('outgroup'):
        s=allele_freqs[n_outgroup,:]
    elif scale_type.startswith('average'):
        s=mean(allele_freqs, axis=0)
    else:
        scaler=1.0
    if scale_type.endswith('product'):
        mu=mean(s)
        scaler=mu*(1.0-mu)
    elif scale_type.endswith('sum'):
        scaler=mean(s*(1.0-s))
    return scaler

def avg_var(ps):
    return sum((p*(1-p) for p in ps))/float(len(ps))

def heterogeneity(allele_frequency, pop_size, type_of_scaling='unbiased'):
    if type_of_scaling=='mle':
        mult=2.0/float(len(allele_frequency))#*float(pop_size)/float(pop_size-1)
    else:
        mult=2.0/float(len(allele_frequency))*float(pop_size)/float(pop_size-1)
    return sum([p*(1-p) for p in allele_frequency])*mult

def B(allele_frequency, pop_size, type_of_scaling='unbiased'):
    return heterogeneity(allele_frequency, pop_size, type_of_scaling)/2.0/float(pop_size)

def adjuster(Bs):
    m=len(Bs)
    res=diag(array(Bs))
    res=res-array(Bs)/m
    res=(res.T-array(Bs)/m).T
    res=res+sum(Bs)/m**2
    return res
    

def bias_correction(m, p, pop_sizes, n_outgroup=None, type_of_scaling='unbiased'):
    #pop sizes are the number of chromosome for each SNP. It is also called the haploid population size
    Bs=[B(prow, pop_size, type_of_scaling=type_of_scaling) for prow, pop_size in zip(p, pop_sizes)]
    print 'Bs',Bs
    adjusting_matrix=adjuster(Bs)
    print 'adjusting matrix',adjusting_matrix
    print 'm',m
    #if n_outgroup is not None:
    #    adjusting_matrix[n_outgroup,:]=0
    #    adjusting_matrix[:,n_outgroup]=0
    print 'adjusting matrix',adjusting_matrix
    res=m-adjusting_matrix
    print 'm-adjusting', res
    from reduce_covariance import reduce_covariance
    print 'mreduced', reduce_covariance(m, n_outgroup)
    print 'adjustingreduced', reduce_covariance(adjusting_matrix, n_outgroup)
    print 'mreduced -adjusting reduced', reduce_covariance(m, n_outgroup)-reduce_covariance(adjusting_matrix, n_outgroup)
    print '(m -adjusting) reduced', reduce_covariance(res, n_outgroup)
    return adjusting_matrix

def other_bias_correction(m,p,pop_sizes,n_outgroup):
    counter=0
    p_vars=[]
    for i in range(len(pop_sizes)):
        if i==n_outgroup:
            p0_var=avg_var(p[n_outgroup,:])/(pop_sizes[n_outgroup])
        else:
            p_vars.append(avg_var(p[i,:])/(pop_sizes[i]))
    adjusting_matrix=diag(array(p_vars))+p0_var
    print 'adjusting_matrix', adjusting_matrix
    print 'm', m
    res=m-adjusting_matrix
    print 'm-adjusting_matrix', res
    return res

class ScaledEstimator(Estimator):
    
    def __init__(self, reducer='',
                 reduce_method=['outgroup','average','None'],
                 scaling=['None','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product','Jade','Jade-o'],
                 full_nodes=None,
                 reduce_also=True,
                 variance_correction=['None','unbiased','mle'],
                 jade_cutoff=1e-5,
                 bias_c_weight='default',
                 names=None):
        super(ScaledEstimator, self).__init__(outgroup_name=reducer, 
                                              full_nodes=full_nodes, 
                                              reduce_also=reduce_also, 
                                              names=names)
        self.scaling=initor(scaling)
        self.variance_correction=initor(variance_correction)
        self.jade_cutoff=jade_cutoff
        self.bias_c_weight=default_scale_dic[bias_c_weight]
        self.reduce_method=reduce_method
        
    def subtract_ancestral_and_get_outgroup(self,p):
        if self.reduce_method=='outgroup':
            n_outgroup=self.get_reduce_index()
            #print n_outgroup
            return p-p[n_outgroup,:], n_outgroup
        elif self.reduce_method=='average':
            n_outgroup=n_outgroup=self.get_reduce_index()
            total_mean2=mean(p, axis=0)
            return p2-total_mean2, n_outgroup
        else:
            return p, None
        
    def __call__(self, xs, ns):
        p=xs/ns
        return self.estimate_from_p(xs,ns)
        
    def estimate_from_p(self, p):
        p=reorder_rows(p, self.names, self.nodes)
        
        p2,n_outgroup = self.subtract_ancestral_and_get_outgroup(p)
        
        
    
        if self.scaling=='Jade':
            mu=mean(p, axis=0)
            
            i=array([v > self.jade_cutoff and v<1.0-self.jade_cutoff for v in mu ])
            p2=p2[:,i]/sqrt(mu[i]*(1.0-mu[i]))
        elif method_of_weighing_alleles=='Jade-o':
            mu=p[n_outgroup,:]
            
            i=array([v > self.jade_cutoff and v<1.0-self.jade_cutoff for v in mu ])
            #p=p[:,i]
            p2=p2[:,i]/sqrt(mu[i]*(1.0-mu[i]))
        
        m=p2.dot(p2.T)/p2.shape[1]
        m=m/m_scaler(m, self.scaling, p, n_outgroup)
        
    
        if reduce_also:
            m=reduce_covariance(m, n_outgroup)
            if self.variance_correction!='None':
                warnings.warn('assuming the same population size for all SNPs', UserWarning)
                pop_sizes=mean(ns, axis=1)
                m=other_bias_correction(m,p, pop_sizes,n_outgroup)
        elif self.variance_correction!='None':
            warnings.warn('assuming the same population size for all SNPs', UserWarning)
            pop_sizes=mean(ns, axis=1)
            changer=bias_correction(m,p, pop_sizes,n_outgroup, type_of_scaling=self.variance_correction)/self.bias_c_weight(m, self.bias_c_weight, p, n_outgroup)
            m=m/m_scaler(m, method_of_weighing_alleles, p, n_outgroup)
            #print 'm',reduce_covariance(m,n_outgroup)
            #print 'changer', reduce_covariance(changer, n_outgroup)
            m-=changer
        
        return m