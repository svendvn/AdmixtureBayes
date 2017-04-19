from Rproposal_admix import addadmix_class, deladmix_class
from Rproposal_regraft import regraft_class
from Rproposal_rescale import rescale_class
from numpy.random import choice
from Rtree_operations import get_number_of_admixes

class new_node_naming_policy(object):
    
    def __init__(self, n=0):
        self.n=0
        
    def next_nodes(self, no_nodes):
        #print self.n
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
        self.props=[addadmix_class(), deladmix_class(), regraft_class(), rescale_class()]
        self.params=[None, None, None, [0.1]]
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        if get_number_of_admixes(tree)==0 and index<=1:
            index=0
            backj=0.25
            forwj=0.5
        elif get_number_of_admixes(tree)==1 and index==1:
            backj=0.5
            forwj=0.25
        else:
            backj=0.25
            forwj=0.25
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,forwj,backj

    def adapt(self,mhr):
        pass
    
    def extract_new_values(self):
        information={}
        information['n']=self.node_naming.n
        return information
    
    def wear_new_values(self, information):
        self.node_naming.n=information['n']
    
class no_admix_proposal(object):
    
    def __init__(self):
        self.props=[regraft_class(), rescale_class()]
        self.params=[None, [0.5]]
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
        return new_tree,forward,backward,1,0.5,0.5
    
    
    
    
    
    
if __name__=='__main__':
    bmp=basic_meta_proposal()
    from Rcatalogue_of_trees import *
    from tree_plotting import pretty_string
    tree=tree_good
    for _ in xrange(50):
        tree=bmp(tree_good)[0]
    print pretty_string(tree)