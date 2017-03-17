from Rproposal_admix import addadmix_class, deladmix_class
from Rproposal_regraft import regraft_class
from Rproposal_rescale import rescale_class
from numpy.random import choice

class new_node_naming_policy(object):
    
    def __init__(self, n=0):
        self.n=0
        
    def next_nodes(self, no_nodes):
        if no_nodes==2:
            return ['x'+str(self.n)+a for a in ['a','b']]
            self.n+=1    
        elif no_nodes==1:
            return 'x'+str(self.n)
            self.n+=1    
        else:
            return ''

class basic_meta_proposal(object):
    
    def __init__(self):
        self.props=[addadmix_class(),deladmix_class(), regraft_class(), rescale_class()]
        self.params=[None, None, None, 0.01]
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        if self.params[index] is not None:
            new_tree, forward, backward =self.props[index](tree, self.params[index], pks=pks)
        else:
            new_tree, forward, backward =self.props[index](tree, pks=pks)
        return new_tree,forward,backward,1,0.25,0.25
    
    
    
    
    
if __name__=='__main__':
    bmp=basic_meta_proposal()
    from Rcatalogue_of_trees import *
    print bmp(tree_good)