import subprocess
from numpy import ix_, array
from reduce_covariance import reduce_covariance

def split_around_delimiter(text, delimiter=','):
    #print text
    return map(int,text.split(delimiter))

def get_muhat(filename):
    minor_alleles=0
    major_alleles=0
    with open(filename, 'r') as f:
        f.readline() #skipping first line with population names
        for fl in f.readlines():
            for pair in fl.split():
                adds=split_around_delimiter(pair)
                minor_alleles+=adds[0]
                major_alleles+=adds[1]
    muhat=float(minor_alleles*major_alleles)/float(major_alleles+minor_alleles)**2
    print 'muhat2', muhat
    return muhat
        
        

def read_data(filename, outgroup='', blocksize=1, nodes=None, noss=False, normalize=False, reduce_also=False, reducer='', return_muhat=False, outfile='tmp'):
    
    args=['treemix', '-i', filename, '-o', outfile,  '-m', '0','-k', str(blocksize)]
    if noss:
        args.append('-noss')
    print args
    subprocess.call(args)#, shell=True)
    args2=['cp', outfile+'.cov', outfile+'.cov.tmp',';','gunzip', '-f',outfile+'.cov.gz', ';', 'mv', outfile+'.cov.tmp', outfile+'.cov']
    args2=['gunzip', '-f',outfile+'.cov.gz']
    
    print 'args2', args2
    subprocess.call(args2)
    
    with open(outfile+'.cov', 'r') as f:
        cats=f.readline().split()
        if nodes==None:
            nodes=cats
        res=[]
        for l in f.readlines():
            numbers=map(float,l.split()[1:])
            res.append(numbers)
            
    mapping={val:key for key, val in enumerate(cats)}
    print nodes
    print mapping
    new_order=[mapping[node] for node in nodes]
    res=array(res)
    if normalize:
        #uncompressing
        uncompressed_file=filename[:-3]
        args3=['gunzip', '-f',filename]
        subprocess.call(args3)
        muhat=get_muhat(uncompressed_file)
        res=res/muhat #remo
    res=res[:,new_order][new_order]
    
    if reduce_also:
        if isinstance(reducer, basestring):
            reduce_index=mapping[reducer]
        else:
            reduce_index=reducer
        res=reduce_covariance(res, reduce_index)
    
    if return_muhat:   
        return res, muhat
    return res

if __name__=='__main__':
    
    print read_data('example1.treemix_in.gz')
    
        
    
    
    
    