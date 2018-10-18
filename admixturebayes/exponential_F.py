from scipy.stats import f
from math import exp


def rvs(t,delta_L, df1=100, df2=100):
    random_part= f.rvs(df1,df2)
    add=delta_L+t
    non_random_part=0
    if add>1:
        non_random_part=add
    else:
        non_random_part=exp(add-1)
    print random_part, 1/random_part
    return t+random_part*(non_random_part-t)

def logpdf(x,t,delta_L, df1=100, df2=100):
    add=delta_L+t
    non_random_part=0
    if add>1:
        non_random_part=add
    else:
        non_random_part=exp(add-1)
    random_part = (x-t)/(non_random_part-t)
    print random_part, 1/random_part
    return f.logpdf(random_part, df1,df2)

if __name__=='__main__':
    delta_L=1.5
    t=0.2
    res=[]
    impossibles=0
    for _ in xrange(100):
        x= rvs(t,delta_L)
        lpdf= logpdf(x,t,delta_L)
        lpdf_back= logpdf(t,x,-delta_L)
        diff=-lpdf+lpdf_back
        res.append(diff)
        
    from numpy import mean
    print min(res), max(res)
    print mean([min(1,exp(l)) for l in res if l<5])
    print len([l for l in res if l>5])
    
    
    