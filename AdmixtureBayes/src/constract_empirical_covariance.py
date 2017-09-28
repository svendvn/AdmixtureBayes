import subprocess
from numpy import array, mean, zeros
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

def treemix_to_cov(filename='treemix_in.txt.gz', 
                   reduce_method=['no', 'average', 'outgroup'], 
                   reducer='', 
                   noss='not implemented', 
                   blocksize='unused', 
                   nodes=None,
                   method_of_weighing_alleles=['None', 'Jade']):
#unzip
    new_filename=make_uncompressed_copy(filename)
    
    allele_counts, names, pop_sizes, minors, total_sum= read_freqs(new_filename)
    
    p=array(allele_counts)
    p=p.T
    p=reorder_rows(p, names, nodes)
    
    if not isinstance(reduce_method, basestring):
        reduce_method=reduce_method[0]
    
    if reduce_method=='outgroup':
        n_outgroup=next((n for n, e in enumerate(names) if e==reducer))
        print n_outgroup
        p=p-p[n_outgroup,:]
    elif reduce_method=='average':
        n_outgroup=next((n for n, e in enumerate(names) if e==reducer))
        print n_outgroup
        p=p-mean(p, axis=0)
    elif reduce_method:
        pass
    
    
        
    m=p.dot(p.T)/p.shape[1]
    
    if noss:
        m=bias_correction(m,p, pop_sizes)
    #print m
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
    res=zeros((m,m))
    res=res-array(Bs)/m
    res=res-array(Bs).T/m
    
    

def bias_correction(m, p, pop_sizes):
    #pop sizes are the number of chromosome for each SNP. It is also called the haploid population size
    Bs=[B(prow, pop_size) for prow, pop_size in zip(p, pop_sizes)]
    adjusting_matrix=adjuster(Bs)
    
    