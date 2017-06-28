import subprocess
from numpy import ix_, array

def read_data(filename, outgroup='Yoruba', blocksize=1, nodes=None, noss=False):
    
    args=['treemix', '-i', filename, '-o', 'tmp', '-root', outgroup, '-m', '0','-k', str(blocksize)]
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
    res=res[:, new_order][new_order]
    return res

if __name__=='__main__':
    
    print read_data('example1.treemix_in.gz')
    
        
    
    
    
    