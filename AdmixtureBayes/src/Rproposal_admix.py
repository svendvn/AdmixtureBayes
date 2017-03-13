from numpy.random import choice, random, exponential
from copy import deepcopy
from scipy.special import binom
from Rtree_operations import get_number_of_admixes, node_is_admixture, insert_admixture_node_halfly, get_descendants_and_rest, graft
from random import getrandbits


def _get_possible_starters(tree):
    res=[]
    for key, node in tree.items():
        if node_is_admixture(node):
            res.extend([(key, 0),(key,1)])
        else:
            res.append((key,0))

def _get_possible_sources(tree, children, other):
    '''
    returns the keys of all non-rooted nodes, that can be grafted into.
    '''
    res=[]
    for oth in other:
        if oth in tree:
            res.append((oth,0))
            if tree[oth][1] is not None:
                res.append((oth,1))
    for child in children:
        if child in tree:
            if tree[child][0] in other:
                res.append((child,0))
            elif tree[child][1] is not None and tree[child][1] in other:
                res.append((child,1))
    return res


def addadmix(tree,pks={}):
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
    
    possible_nodes=_get_possible_starters(tree)
        
    new_tree= deepcopy(tree)
    sink_key, sink_branch=choice(possible_nodes, 1)[0]
    children, other= get_descendants_and_rest(tree, sink_key)
    candidates=_get_possible_sources(new_tree, children, other)+[('r',0)]
    ch= choice(len(candidates),1)[0]
    source_key, source_branch=candidates[ch]
    print 'sink', sink_key
    print 'source', (source_key,source_branch)
    print 'new_tree',new_tree
    new_tree, forward_backward= insert_admix(new_tree, regrafter, recipient_key, new_node=new_node, which_branch=recipient_branch)
    _, new_other =  get_descendants_and_rest(new_tree, regrafter)
    
def get_admixture_branch_length(x=None):
    if x is None:
        return 0,1
    else:
        return 1
    

def insert_admix(tree, source_key, source_branch, sink_key, sink_branch, source_name=None, sink_name=None):
    if source_key=='r':
        u1=exponential()
    else:
        u1=random()
    u2=random()
    t1,q1=get_admixture_branch_length()
    if sink_name is None:
        sink_name=str(getrandbits(8)).strip()
    tree=insert_admixture_node_halfly(tree, source_key, source_branch, u2, admix_b_length=t1, new_node=sink_name)
    tree=graft(tree, regrafter, recipient_key, new_node=new_node, which_branch=recipient_branch)

def deladmix(tree,pks={}):
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
    