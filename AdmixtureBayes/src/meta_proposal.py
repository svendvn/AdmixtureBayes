from Rproposal_admix import addadmix_class, deladmix_class
from Rproposal_regraft import regraft_class
from Rproposal_rescale import rescale_class
from numpy.random import choice

class new_node_naming_policy(object):
    
    def __init__(self, n=0):
        self.n=0
        
    def next_nodes(self, no_nodes):
        if no_nodes==2:
            self.n+=1 
            return ['x'+str(self.n)+a for a in ['a','b']]
        elif no_nodes==1:
            self.n+=1  
            return 'x'+str(self.n)
              
        else:
            return ''

class basic_meta_proposal(object):
    
    def __init__(self):
        self.props=[addadmix_class(),deladmix_class(), regraft_class(), rescale_class()]
        self.params=[None, None, None, [0.01]]
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,0.25,0.25
    
    
    
    
    
if __name__=='__main__':
    bmp=basic_meta_proposal()
    from Rcatalogue_of_trees import *
    from tree_plotting import pretty_string
    tree=tree_good
    for _ in xrange(50):
        tree=bmp(tree_good)[0]
    print pretty_string(tree)