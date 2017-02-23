from numpy.random import choice
from proposal_admix import addadmix, deladmix
from proposal_rescale import rescale

all_props=[addadmix,deladmix,rescale]
all_props_names=['addadmix','deladmix','rescale']

def prop_flat(x,type_chooser, pks={}):
    '''
    Input x: tree
    output tuple:
        output  z: new tree
        output g1: density of jumping forwards
        output g2: density of jumping backwards
        output Jh: the Jacobian of the transition(necessary for Reversible Jump MCMC).
        output j1: probability of trying forward jump
        output j2: probability of trying backward jump
    '''
    
    jump_type=choice(len(all_props),1)[0]
    j1=j2=1.0/float(len(all_props))
    pks['proposal_type']=all_props_names[jump_type]
    pks['j1_type']=j1
    pks['j2_type']=j2
    
    return (all_props[jump_type](x,pks)+(j1,j2))
    
    