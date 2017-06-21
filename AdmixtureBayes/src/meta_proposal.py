from Rproposal_admix import addadmix_class, deladmix_class
from Rproposal_regraft import regraft_class
from Rproposal_rescale import rescale_class
from Rproposal_sliding_regraft import sliding_regraft_class, sliding_regraft_class_resimulate
from Rproposal_rescale_marginally import rescale_marginally_class
from Rproposal_sliding_rescale import sliding_rescale_class
from numpy.random import choice
from Rtree_operations import get_number_of_admixes
from math import exp

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
        self.props=[addadmix_class(), deladmix_class(), regraft_class(), rescale_class(), rescale_marginally_class()]
        self.params=[None, None, None, [0.1], [0.1]]
        self.node_naming=new_node_naming_policy()
        
    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        if get_number_of_admixes(tree)==0 and index<=1:
            index=0
            backj=0.2
            forwj=0.4
        elif get_number_of_admixes(tree)==1 and index==1:
            backj=0.4
            forwj=0.2
        else:
            backj=0.2
            forwj=0.2
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,forwj,backj

    def adapt(self,mhr, u, post_new, post, temperature):
        pass
    
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information     
    
    def wear_exportable_state(self, information):
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
    
    def adapt(self,mhr, u, post_new, post, temperature):
        pass
    
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']

class adaptive_proposal_no_admix(object):
    
    def __init__(self):
        self.props=[sliding_regraft_class(), rescale_class()]
        start_value_of_sigma=0.1
        start_value_of_slider=0.1
        self.node_naming=new_node_naming_policy()
        self.recently_called_type=None
        self.regraft_count=10
        self.rescale_count=10
        self.multiplier=10
        self.desired_mhr=0.234
        self.alpha=0.9
        self.params=[[start_value_of_slider], [start_value_of_sigma]]

    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        backj=0.5
        forwj=0.5
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        self.recently_called_type=self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,forwj,backj
    
    def adapt(self,mhr, u, post_new, post, temperature):
        if self.recently_called_type=='rescale':
            new_val, self.rescale_count= standard_update(self.rescale_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[1][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='rescale')
            self.params[1]=[new_val]
        if self.recently_called_type=='sliding_regraft':
            new_val, self.regraft_count= standard_update(self.regraft_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[0][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='regraft_slider',
                                                         max_val=15.0)
            self.params[0]=[new_val]
            
            
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        #information['params']=self.params
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']
        #self.params=information['params']

class adaptive_proposal(object):
    
    def __init__(self, resimulate_regrafted_branch_length=False):
        self.props=[addadmix_class(), deladmix_class(), sliding_regraft_class(), rescale_class(), sliding_rescale_class()]
        if resimulate_regrafted_branch_length:
            self.props[2]=sliding_regraft_class_resimulate(resimulate_regrafted_branch_length)
        start_value_of_sigma=0.1
        start_value_of_slider=0.1
        start_value_of_marginal_sigma=0.1
        self.node_naming=new_node_naming_policy()
        self.recently_called_type=None
        self.regraft_count=10
        self.rescale_count=10
        self.multiplier=10
        self.desired_mhr=0.234
        self.alpha=0.9
        self.params=[None, None, [start_value_of_slider], [start_value_of_sigma], [start_value_of_marginal_sigma]]
        
    def __call__(self, tree, pks={}):
        index=choice(len(self.props),1)[0]
        if get_number_of_admixes(tree)==0 and index<=1:
            index=0
            backj=1.0/len(self.props)
            forwj=2.0/len(self.props)
        elif get_number_of_admixes(tree)==1 and index==1:
            backj=2.0/len(self.props)
            forwj=1.0/len(self.props)
        else:
            backj=1.0/len(self.props)
            forwj=1.0/len(self.props)
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        pks['proposal_type']= self.props[index].proposal_name
        self.recently_called_type=self.props[index].proposal_name
        args=[]
        if names:
            args.append(names)
        if self.params[index] is not None:
            args.extend(self.params[index])
        new_tree, forward, backward =self.props[index](tree, *args, pks=pks)
        return new_tree,forward,backward,1,forwj,backj

    def adapt(self, mhr, u, post_new, post, temperature):
        if self.recently_called_type=='rescale':
            new_val, self.rescale_count= standard_update(self.rescale_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[3][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='rescale')
            self.params[3]=[new_val]
        if self.recently_called_type=='sliding_regraft':
            new_val, self.regraft_count= standard_update(self.regraft_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[2][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='regraft_slider',
                                                         max_val=15.0)
            self.params[2]=[new_val]
            
        if self.recently_called_type=='sliding_rescale':
            new_val, self.regraft_count= standard_update(self.regraft_count, 
                                                         self.multiplier, 
                                                         self.alpha, 
                                                         self.params[4][0], 
                                                         mhr, 
                                                         desired_mhr=self.desired_mhr, 
                                                         verbose=False,
                                                         name='sliding_rescale',
                                                         max_val=15.0)
            self.params[4]=[new_val]
            
            
    
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        #information['params']=self.params
        return information     
    
    def wear_exportable_state(self, information):
        self.node_naming.n=information['n']
        #self.params=information['params']
        
        
def standard_update(count, multiplier, alpha, old_value, mhr, desired_mhr=0.234, verbose=False, max_val=float('inf'), name='value'):
    count+=1
    gamma=multiplier/count**alpha
    change=exp(gamma*(min(1.0,mhr)-desired_mhr))
    value=old_value*change
    value=min(value, max_val)
    if verbose:
        print 'old_'+name+'=',old_value
        print 'mhr=',mhr
        print 'count=',count
        print 'gamma=', gamma
        print 'multiplier=', change
        print 'new_'+name+'=', value
    return value,count
    

    
    
    
    
if __name__=='__main__':
    bmp=basic_meta_proposal()
    from Rcatalogue_of_trees import *
    from tree_plotting import pretty_string
    tree=tree_good
    for _ in xrange(50):
        tree=bmp(tree)[0]
    print pretty_string(tree)