import subprocess
from numpy import ix_, array

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
        
        

def read_data(filename, outgroup='', blocksize=1, nodes=None, noss=False, normalize=True):
    
    if outgroup:
        args=['treemix', '-i', filename, '-o', 'tmp', '-root', outgroup, '-m', '0','-k', str(blocksize)]
    else:
        args=['treemix', '-i', filename, '-o', 'tmp',  '-m', '0','-k', str(blocksize)]
    if noss:
        args.append('-noss')
    print args
    subprocess.call(args)#, shell=True)
    args2=['gunzip', '-k', '-f','tmp.cov.gz']
    subprocess.call(args2)
    
    with open('tmp.cov', 'r') as f:
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
        args3=['gunzip', '-k', '-f',filename]
        subprocess.call(args3)
        res=res/get_muhat(uncompressed_file) #remo
    res=res[:,new_order][new_order]
    return res

if __name__=='__main__':
    
    print read_data('example1.treemix_in.gz')
    
        
    
    
    
    