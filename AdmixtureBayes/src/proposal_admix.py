from numpy.random import random, choice
from copy import deepcopy
from scipy.special import binom
from tree_operations import get_number_of_admixes, get_all_branch_lengths
from random import getrandbits


def _update_branch(b,other,direction,prop,identifier):
    '''
    This function takes a branch, B, of the form [root, destination, branch length 1, admixture event 1, branch length 2, admixture event 2, ..., branch length n]
    and induces an extra admixture event at the relative position, pos_de_new, which are generated in this function. 
    The direction is in DIRECTION and the source/sink is encoded in OTHER. PROP is the admixture proportion.
    '''
    #the times of events on the branch
    times=b[2::2]
    b_length=sum(times)
    u=random()
    pos_de_new=u*b_length
    cumtime=0
    #print pos_de_new,b_length
    for n,time in enumerate(times):
        old_cumtime=cumtime
        cumtime+=time
        if cumtime>(pos_de_new-1e-7):
            insertion=[pos_de_new-old_cumtime, [direction, other, prop, identifier], cumtime-pos_de_new]
            break
    return b[:(2+n*2)] + insertion + b[(2+n*2+1):], b_length, u

def _jacobian(c1,c2,u1,u2,w):
    '''
    This function returns the jacobian of h(c1,c2,u1,u2,w)=(c1*u1, c1*(1-u1), c2*u2, c2*(1-u2), 0.5*w)
    '''
    return -0.5*c1*c2

def addadmix(tree):
    '''
    This proposal adds an admixture to the tree. There are a lot of free parameters but only 5 are in play here:
        c1: the branch length of the source population
        c2: the branch legnth of the population genes are migrating into. 
        u1: the position of the admixture source on the source population branch
        u2: the position of the admixture destination on the sink population branch
        w: the admixture proportion.
        The connecting function, h, (see Green & Hastie 2009) is
            h(c1,c2,u1,u2,w)=(c1*u1, c1*(1-u1), c2*u2, c2*(1-u2), 0.5*w)
    '''
    
    no_admixs=get_number_of_admixes(tree)
    
    i1,i2=choice(len(tree), 2, replace=False) #taking two random admixtures.
    
    cop=deepcopy(tree)     
    
    identifier=getrandbits(128)
    
    cop[i1],c1,u1=_update_branch(cop[i1], cop[i2][1],">",None,identifier)
    w=random() #admixture proportion.
    cop[i2],c2,u2=_update_branch(cop[i2], cop[i1][1],"<", w,identifier)
    
    absolut_jacobian=1#abs(_jacobian(c1,c2,u1,u2,w))
        
    return cop, 1,1,1#1.0/(binom(len(cop),2)*2), 1.0/float(no_admixs+1),absolut_jacobian 


def deladmix(tree):
    '''
    Reversible Jump MCMC transition which removes a random admixture branch. This is the reverse of the other proposal in this file. 
    '''
    
    #necessary for calculation of transition probabilities.
    no_admixs=get_number_of_admixes(tree)
    
    #if there are no admixture events there is nothing to erase and we should just move on
    if no_admixs==0:
        return tree,1,1,1
    
    #making copy that we can erase branches from. 
    cop=deepcopy(tree)
    
    #get a random admixture event
    i=choice(no_admixs*2, 1)[0]
    for n,branch in enumerate(cop):
        no_admixes_in_branch=(len(branch)-3)/2
        if no_admixes_in_branch>i:
            source=branch[1]
            admixture_partner_branch_key=branch[3+i*2][1]
            admixture_partner_identifier=branch[3+i*2][3]
            c1=branch[(3+i*2-1)]+branch[(3+i*2+1)]
            cop[n]=branch[:(3+i*2-1)]+[c1]+branch[(3+i*2+2):] #removing the branch
            break
        else:
            i-=no_admixes_in_branch
    for n,branch in enumerate(cop):
        if admixture_partner_branch_key==branch[1]:
            no_admixes_in_branch=(len(branch)-3)/2
            for j in range(0,no_admixes_in_branch):
                if branch[3+j*2][1]==source and branch[3+2*j][3]==admixture_partner_identifier:
                    c2=branch[(3+j*2-1)]+branch[(3+j*2+1)]
                    cop[n]=branch[:(3+j*2-1)]+[c2]+branch[(3+j*2+2):] #removing the receiver branch
                    absolute_jacobian=1.0#/abs(_jacobian(c1,c2,0,0,0)) #u1,u2 and w are not extracted. 
                    return cop,1,1,1#1.0/float(no_admixs), 1.0/(binom(len(cop),2)*2), absolute_jacobian


class Tester():
    
    def __init__(self, tree=None):
        if tree is None:
            tree=[["r","s1",0.3,[">","s3",None,32123], 0.1], 
                   ["s2s3","s2",0.2],
                   ["s2s3","s3",0.15, ["<","s1",0.44,32123],0.05],
                   ["r","s2s3",0.2]]
        self.tree=tree
        
    def test_proposals(self):
        a=addadmix(self.tree)[0]
        if get_number_of_admixes(a)==2:
            print "ADDED ADMIXTURE: TRUE"
        else:
            print "ADDED ADMIXTURE: FALSE"
        b=deladmix(a)[0]
        c=deladmix(b)[0]
        if get_number_of_admixes(b)==1 and get_number_of_admixes(c)==0:
            print "DELETED ADMIXTURE: TRUE"
        else:
            print "DELETED ADMIXTURE: FALSE"
        d=deladmix(c)[0]
        if get_number_of_admixes(d)==0:
            print "HANDLED DELETION OF TRIVIAL ADMIXTURE: TRUE"
        else:
            print "HANDLED DELETION OF TRIVIAL ADMIXTURE: FALSE"

    def test_maintenance_of_branch_length(self):
        before=get_all_branch_lengths(self.tree)
        a=addadmix(self.tree)[0]
        after=get_all_branch_lengths(a)
        
        
        
        if set(before.keys()) == set(after.keys()):
            print "ADD ADMIXTURE MAINTAINS BRANCH TOPOLOGY: TRUE"
        else:
            print "ADD ADMIXTURE MAINTAINS BRANCH TOPOLOGY: FALSE"
            return
        
        for branch,length in before.items():
            if abs(length-after[branch])>1e-7:
                print "ADD ADMIXTURE MAINTAINS BRANCH LENGTHS: FALSE"
                return
        print "ADD ADMIXTURE MAINTAINS BRANCH LENGTHS: TRUE"
        
        before=after
        b=deladmix(a)[0]
        after=get_all_branch_lengths(b)
        
        if set(before.keys()) == set(after.keys()):
            print "DELETE ADMIXTURE MAINTAINS BRANCH TOPOLOGY: TRUE"
        else:
            print "DELETE ADMIXTURE MAINTAINS BRANCH TOPOLOGY: FALSE"
            return
        
        for branch,length in before.items():
            if abs(length-after[branch])>1e-9:
                print "DELETE ADMIXTURE MAINTAINS BRANCH LENGTHS: FALSE"
                return
        print "DELETE ADMIXTURE MAINTAINS BRANCH LENGTHS: TRUE"
               

if __name__=="__main__":
    t=Tester()
    t.test_proposals()
    t.test_maintenance_of_branch_length()
    