import newick
from copy import deepcopy
from numpy.random import normal, choice, random
from scipy.special import binom
from numpy import random as npr
from numpy import cumsum
from tree_operations import *



def deladmix(tree):
    '''
    Reversible Jump MCMC transition which removes a random admixture branch. 
    '''
    
    
    cop=deepcopy(tree)
    
    #necessary for calculation of transition probabilities.
    tot_length=sum(map(len, cop))
    no_admixs=(tot_length-len(cop)*3)/4
    
    #if there are no admixture events there is nothing to erase and we should just move on
    if no_admixs==0:
        return cop,1,1
    
    #get a random admixture event
    i=choice(no_admixs*2, 1)[0]
    for branch in cop:
        no_admixes_in_branch=(len(branch)-3)/2
        if no_admixes_in_branch>i:
            source=branch[1]
            admixture_partner=branch[3+i*2][1]
            branch=branch[:(3+i*2-1)]+[branch[(3+i*2-1)]+branch[(3+i*2+1)]]+branch[(3+i*2+1):] #removing the branch
            break
        else:
            i-=no_admixes_in_branch
    for branch in cop:
        if admixture_partner==branch[1]:
            no_admixes_in_branch=(len(branch)-3)/2
            for j in range(0,no_admixes_in_branch):
                if branch[3+j*2][1]==source:
                    branch=branch[:(3+j*2-1)]+[branch[(3+j*2-1)]+branch[(3+j*2+1)]]+branch[(3+j*2+1):] #removing the branch
                    return cop,1.0/float(no_admixs), 1.0/(binom(len(cop),2)*2)
                
                
def addadmix(tree):
    '''
    This proposal adds an admixture to the tree
    '''
    tot_length=sum(map(len, tree))
    no_admixs=(tot_length-len(tree)*3)/4
    
    i1,i2=npr.choice(len(tree), 2, replace=False) #taking two random admixtures.
    
    cop=deepcopy(tree)
    
    def update_branch(b,other,direction,prop):
        #the times of events on the branch
        times=b[2::2]
        b_length=sum(times)
        pos_de_new=random()*b_length
        cumtime=0
        for n,time in enumerate(times):
            old_cumtime=cumtime
            cumtime+=time
            if cumtime>pos_de_new:
                insertion=[pos_de_new-old_cumtime, [direction, other, prop], cumtime-pos_de_new]
                break
        return b[:(2+n*2)]+insertion+b[(2+n*2+1):]
        
    
    cop[i1]=update_branch(cop[i1], cop[i2][1],">",None)
    cop[i2]=update_branch(cop[i2], cop[i1][1],"<", random()*0.5)
        
    return cop, 1.0/(binom(len(cop),2)*2), 1.0/float(no_admixs+1)
    
    

print rescale(rescale(tree_flatter_list)[0])

print deladmix(tree_flatter_list)
print addadmix(tree_flatter_list)

#visitor=Edges()
#print visitor.get_edges(tree_string)

admixtures=[("3","1",0.2,0.8,0.3,1)]
receivers=["1"]
senders=["3"]





#print newick.tree.parse_tree(, )
    