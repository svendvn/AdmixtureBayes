import subprocess
from numpy import array, mean, zeros, diag, sum, arcsin, sqrt
from reduce_covariance import reduce_covariance

def reorder_rows(p, names,nodes):
    mapping={val:key for key, val in enumerate(names)}
    if nodes is None:
        nodes=names
    new_order=[mapping[node] for node in nodes]
    p=p[new_order,:]
    return p

def read_freqs(new_filename):
    with open(new_filename, 'r') as f:
        names=f.readline().split()
        allele_counts=[]
        pop_sizes=[]
        minors=[]
        total_sum=0
        for n,r in enumerate(f.readlines()):
            minor_majors=r.split()
            freqs=[]
            for minor_major in minor_majors:
                minor, major= map(float,minor_major.split(','))
                freqs.append(float(minor)/float(major+minor))
                total_sum+=major+minor
                minors.append(minor)
                if n==0:
                    pop_sizes.append(major+minor)
            allele_counts.append(freqs)
    return allele_counts, names, pop_sizes, minors, total_sum

def make_uncompressed_copy(filename):
    take_copy_args=['cp', filename, filename+".tmp"]
    move_back_args=['mv', filename+'.tmp', filename]
    args=['gunzip', '-f', filename]
    new_filename='.'.join(filename.split('.')[:-1])
    subprocess.call(take_copy_args)
    subprocess.call(args)
    subprocess.call(move_back_args)
    return new_filename

def scale_m(m, scale_type, allele_freqs, n_outgroup=None):
    if scale_type.startswith('outgroup'):
        s=allele_freqs[n_outgroup,:]
    else:
        s=mean(allele_freqs, axis=0)
    if scale_type.endswith('product'):
        mu=mean(s)
        scaler=mu*(1.0-mu)
    else:
        scaler=mean(s*(1.0-s))
    return m/scaler
        

def treemix_to_cov(filename='treemix_in.txt.gz', 
                   reduce_method=['no', 'average', 'outgroup'], 
                   reducer='', 
                   variance_correction=False, 
                   nodes=None,
                   arcsin_transform=False,
                   method_of_weighing_alleles=['None', 'Jade','outgroup_sum', 'outgroup_product', 'average_outgroup', 'average_product']
                   ):
#unzip
    new_filename=make_uncompressed_copy(filename)
    
    allele_counts, names, pop_sizes, minors, total_sum= read_freqs(new_filename)
    
    p=array(allele_counts)
    p=p.T
    p=reorder_rows(p, names, nodes)
    
    if not isinstance(reduce_method, basestring):
        reduce_method=reduce_method[0]
    if not isinstance(method_of_weighing_alleles, basestring):
        method_of_weighing_alleles=method_of_weighing_alleles[0]
    
    if arcsin_transform:
        p2=arcsin(sqrt(p))*2
    else:
        p2=p
        
    if reduce_method=='outgroup':
        n_outgroup=next((n for n, e in enumerate(nodes) if e==reducer))
        print n_outgroup
        p2=p2-p2[n_outgroup,:]
    elif reduce_method=='average':
        n_outgroup=next((n for n, e in enumerate(nodes) if e==reducer))
        total_mean2=mean(p, axis=0)
        if arcsin_transform:
            total_mean2=arcsin(sqrt(total_mean2))*2
        p2=p2-total_mean2
    else:
        n_outgroup=None
    
    if method_of_weighing_alleles=='Jade':
        mu=mean(p, axis=0)
        p2=p2/sqrt(mu*(1.0-mu))
        
    m=p2.dot(p2.T)/p2.shape[1]
    
    if variance_correction:
        m=bias_correction(m,p, pop_sizes)
        
    if method_of_weighing_alleles != 'None' and method_of_weighing_alleles!='Jade':
        m=scale_m(m, method_of_weighing_alleles, p, n_outgroup)

    if reduce_method != 'no':
        m=reduce_covariance(m, n_outgroup)
        
    return m

def heterogeneity(allele_frequency, pop_size):
    mult=2.0/float(len(allele_frequency))*float(pop_size)/float(pop_size-1)
    return sum([p*(1-p) for p in allele_frequency])*mult

def B(allele_frequency, pop_size):
    return heterogeneity(allele_frequency, pop_size)/2.0/float(pop_size)

def adjuster(Bs):
    m=len(Bs)
    res=diag(array(Bs))
    res=res-array(Bs)/m
    res=(res.T-array(Bs)/m).T
    res=res+sum(Bs)/m**2
    return res
    

def bias_correction(m, p, pop_sizes):
    #pop sizes are the number of chromosome for each SNP. It is also called the haploid population size
    Bs=[B(prow, pop_size) for prow, pop_size in zip(p, pop_sizes)]
    adjusting_matrix=adjuster(Bs)
    return m-adjusting_matrix

if __name__=='__main__':
    
    m=adjuster([0.1,0.2])
    assert sum(m-array([[0.075,-0.075],[-0.075,0.075]]), axis=None)==0, 'adjusting wrong'
    from numpy import set_printoptions
    set_printoptions(precision=5)
    filename='sletmig/_treemix_in.txt.gz'
    nodes=['s'+str(i) for i in range(1,10)]+['out']
    print treemix_to_cov(filename, reduce_method='outgroup', reducer='out', variance_correction=False, nodes=nodes, arcsin_transform=False, method_of_weighing_alleles='outgroup_product')
    
    from load_data import read_data
    
    print read_data(filename, blocksize=1, nodes=nodes, variance_correction=True, normalize=False, reduce_also=True, reducer='out', return_muhat=False)
    