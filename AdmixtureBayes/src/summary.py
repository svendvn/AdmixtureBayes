
## the recipe of a summary function:
#def summary(old_tree,new_tree,old_pos,new_pos,iteration_number, proposal_object):, but all the arguments are encoded as **kwargs to reduce maintenance
#return something

from Rtree_operations import get_number_of_admixes, get_all_branch_lengths

class Summary(object):
       
    def __init__(self, name, pandable=True):
        self.name=name
        self.pandable=pandable
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return str(output)
    
class s_no_admixes(Summary):
    
    def __init__(self):
        super(s_no_admixes,self).__init__('no_admixes')

    def __call__(self, **kwargs):
        old_tree=kwargs['old_tree']
        return get_number_of_admixes(old_tree)

class s_branch_length(Summary):

    def __init__(self):
        super(s_branch_length,self).__init__('branch_length')

    def __call__(self, **kwargs):
        old_tree=kwargs['old_tree']
        return sum(get_all_branch_lengths(old_tree))
    
class s_variable(Summary):
    
    def __init__(self, variable, pandable=True):
        super(s_variable, self).__init__(variable, pandable)

    def __call__(self, **kwargs):
        return kwargs[self.name]

 