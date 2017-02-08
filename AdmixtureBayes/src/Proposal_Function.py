from numpy.random import choice
from proposal_admix import addadmix, deladmix
from proposal_rescale import rescale

mixture=[addadmix,deladmix,rescale]

def prop(x):
    '''
    Input x: tree
    output tuple:
        output z:  new treedensity of jumping back divided by the density of jumping forward
        output g1: density of jumping forwards
        output g2: density of jumping backwards
        output Jh: the jacobian of the transition(necessary for Reversible Jump MCMC).
        output j1: probability of trying forward jump
        output j2: probability of trying backward jump
    '''
    
    jump_type=choice(mixture,1)
    j1=j2=1.0/float(len(mixture))
    
    return(jump_type[0](x)+(j1,j2))
    
    