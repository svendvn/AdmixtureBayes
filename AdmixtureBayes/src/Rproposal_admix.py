from numpy.random import choice, random, exponential
from copy import deepcopy
from scipy.special import binom
from Rtree_operations import (get_number_of_admixes, node_is_admixture, insert_admixture_node_halfly, 
                              get_descendants_and_rest, graft, remove_admix, node_is_non_admixture,
                              make_consistency_checks, parent_is_spouse, halfbrother_is_uncle,
                              parent_is_sibling, other_branch)
from random import getrandbits




def _get_possible_starters(tree):
    res=[]
    for key, node in tree.items():
        if node_is_admixture(node):
            res.extend([(key, 0),(key,1)])
        else:
            res.append((key,0))
    return res

def _get_possible_sources(tree, children, other, sink_key, sink_branch):
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
            if tree[child][0] in other and (child!=sink_key or 0!=sink_branch):
                res.append((child,0))
            elif tree[child][1] is not None and tree[child][1] in other and (child!=sink_key or 1!=sink_branch):
                res.append((child,1))
    return res

class addadmix_class(object):
    
    new_nodes=2
    
    def __call__(self,*args, **kwargs):
        return addadmix(*args, **kwargs)
    
class deladmix_class(object):
    
    new_nodes=0
    
    def __call__(self,*args, **kwargs):
        return deladmix(*args, **kwargs)

def addadmix(tree,new_node_names=None,pks={}):
    '''
    This proposal adds an admixture to the tree. There are a lot of free parameters but only 5 are in play here:
        c1: the branch length of the source population
        c2: the branch length of the population genes are migrating into. 
        u1: the position of the admixture source on the source population branch
        u2: the position of the admixture destination on the sink population branch
        w: the admixture proportion.
        The connecting function, h, (see Green & Hastie 2009) is
            h(c1,c2,u1,u2,w)=(c1*u1, c1*(1-u1), c2*u2, c2*(1-u2), 0.5*w)
    '''
    
    possible_nodes=_get_possible_starters(tree)
        
    new_tree= deepcopy(tree)
    #print possible_nodes
    sink_key, sink_branch=possible_nodes[choice(len(possible_nodes), 1)[0]]
    children, other= get_descendants_and_rest(tree, sink_key)
    candidates=_get_possible_sources(new_tree, children, other, sink_key,sink_branch)+[('r',0)]
    ch= choice(len(candidates),1)[0]
    source_key, source_branch=candidates[ch]
    #print 'children', children
    #print 'candidates', candidates
    #print 'sink', (sink_key, sink_branch)
    #print 'source', (source_key,source_branch)
    #print 'new_tree',new_tree
    if new_node_names is None:
        new_tree, forward_backward= insert_admix(new_tree, source_key, source_branch, sink_key, sink_branch)
    else:
        new_tree, forward_backward= insert_admix(new_tree, source_key, source_branch, sink_key, sink_branch, new_node_names[0], new_node_names[1])
    return new_tree,1,1
    
def get_admixture_branch_length(x=None):
    if x is None:
        return 0.1,1
    else:
        return 1
    
def get_admixture_proportion(x=None):
    if x is None:
        return 0.48,1
    else: 
        return 1
    

def insert_admix(tree, source_key, source_branch, sink_key, sink_branch, source_name=None, sink_name=None):
    if source_key=='r':
        u1=exponential()
    else:
        u1=random()
    u2=random()
    t1,q1=get_admixture_branch_length()
    u3,q3=get_admixture_proportion()
    if sink_name is None:
        sink_name=str(getrandbits(8)).strip()
    if source_name is None:
        source_name=str(getrandbits(8)).strip()
    tree=insert_admixture_node_halfly(tree, sink_key, sink_branch, u2, admix_b_length=t1, new_node_name=sink_name, admixture_proportion= u3)
    #print 'tree after inserting admixture', tree
    tree=graft(tree, sink_name, source_key, u1, source_name, source_branch, remove_branch=1)
    #print 'tree after grafting', tree
    return tree,None


def deladmix(tree,pks={}):
    '''
    Reversible Jump MCMC transition which removes a random admixture branch. This is the reverse of the other proposal in this file. 
    '''
    
    #necessary for calculation of transition probabilities.
    no_admixs=get_number_of_admixes(tree)
    
    #making copy that we can erase branches from. 
    cop=deepcopy(tree)
    
    candidates=_get_removable_admixture_branches(cop)
    if len(candidates)==0:
        return tree
    remove_key, remove_branch = candidates[choice(len(candidates),1)[0]]
    print 'remove', (remove_key, remove_branch)
    
    return remove_admix(tree, remove_key, remove_branch)[0],1,1

def _check_node(tree,key,direction):
    parent_key=tree[key][direction]
    return ((parent_key=='r' or node_is_non_admixture(tree[parent_key])) and 
            not parent_is_spouse(tree, key, other_branch(direction)) and
            (parent_key=='r' or not halfbrother_is_uncle(tree, key, parent_key)) and
            (parent_key=='r' or not (parent_is_spouse(tree,key,direction) and parent_is_sibling(tree, key, direction))))
        
def _get_removable_admixture_branches(tree):
    res=[]
    for key, node in tree.items():
        if node_is_admixture(node):
            if _check_node(tree, key, 0):
                res.append((key,0))
            if _check_node(tree, key, 1):
                res.append((key, 1))
    return res


class Tester():
    
    def __init__(self, tree):
        self.tree=tree
        self.no_admixes=get_number_of_admixes(self.tree)
    
    def many_admixes(self, n=100):
        for i in xrange(n):
            self.tree=addadmix(self.tree, new_node_names=['n'+str(i)+a for a in ['o','n']])
            #plot_graph(self.tree)
            if self.no_admixes+1==get_number_of_admixes(self.tree):
                print 'INCREASED NUMBER OF ADMIXTURES BY ONE= '+'TRUE'
            else:
                print 'INCREASED NUMBER OF ADMIXTURES BY ONE= '+'FALSE'
            self.no_admixes=get_number_of_admixes(self.tree)
            ad=make_consistency_checks(self.tree, ['s1','s2','s3','s4'])
            if not ad[0]:
                print ad
                plot_graph(self.tree, drawing_name='bad.png')
                break
        plot_as_directed_graph(self.tree)
            
    def alternate_admixes(self, n=1000):
        for i in xrange(n):
            old_tree=deepcopy(self.tree)
            self.tree=addadmix(self.tree, new_node_names=['n'+str(i)+a for a in ['o','n']])
            if self.no_admixes+1==get_number_of_admixes(self.tree):
                print 'INCREASED NUMBER OF ADMIXTURES BY ONE= '+'TRUE'
            else:
                print 'INCREASED NUMBER OF ADMIXTURES BY ONE= '+'FALSE'
            ad=make_consistency_checks(self.tree, ['s1','s2','s3','s4'])
            if not ad[0]:
                print ad
                plot_graph(old_tree, drawing_name='good.png')
                plot_graph(self.tree, drawing_name='bad.png')
                break
            #plot_graph(self.tree)
            old_tree=deepcopy(self.tree)
            self.tree=deladmix(self.tree)
            if self.no_admixes==get_number_of_admixes(self.tree):
                print 'DECREASED NUMBER OF ADMIXTURES BY ONE= '+'TRUE'
            else:
                print 'DECREASED NUMBER OF ADMIXTURES BY ONE= '+'FALSE'
            #plot_graph(self.tree)
            print self.tree
            ad=make_consistency_checks(self.tree, ['s1','s2','s3','s4'])
            if not ad[0]:
                print ad
                plot_graph(old_tree, drawing_name='good.png')
                plot_graph(self.tree, drawing_name='bad.png')
                deladmix(old_tree)
                break
    

if __name__=="__main__":
    from tree_plotting import plot_as_directed_graph, plot_graph
    import Rtree_operations
    #plot_graph(Rtree_operations.tree_on_the_border2_with_children)
    t=Tester(Rtree_operations.tree_on_the_border2_with_children)
    t.many_admixes(10)
    