
## the recipe of a summary function:
#def summary(old_tree,new_tree,old_pos,new_pos,iteration_number, proposal_object):, but all the arguments are encoded as **kwargs to reduce maintenance
#return something

from tree_operations import get_number_of_admixes, get_all_branch_lengths

class Summary():
    
    name='default'
    
    def __call__(self, **kwargs):
        pass
    
    def string_answer(self, output):
        return str(output)
    

class s_posterior(Summary):

    name='posterior'

    def __call__(self, **kwargs):
        return kwargs['old_pos']
    
class s_no_admixes(Summary):

    name='no_admixes'

    def __call__(self, **kwargs):
        old_tree=kwargs['old_tree']
        return get_number_of_admixes(old_tree)

class s_branch_length(Summary):

    name='total_length'

    def __call__(self, **kwargs):
        old_tree=kwargs['old_tree']
        return sum(get_all_branch_lengths(old_tree))


print False and (8%0)
 