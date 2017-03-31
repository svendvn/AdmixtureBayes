import matplotlib.pyplot as plt
from Rtree_operations import get_number_of_admixes, get_all_branch_lengths


class Summary(object):
       
    def __init__(self, name, pandable=True):
        self.name=name
        self.pandable=pandable
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return str(output)
    
    def make_trajectory(self, data, **kwargs):
        return plt.plot(data[self.name],**kwargs)
    
    def make_histogram(self, data, **kwargs):
        return plt.hist(data[self.name], fc=(1, 0, 0, 0.5), normed=True, **kwargs)
    
    def add_histogram(self, a, **kwargs):
        return plt.hist(a,fc=(0, 1, 0, 0.5), normed=True, **kwargs)
    
class s_no_admixes(Summary):
    
    def __init__(self):
        super(s_no_admixes,self).__init__('no_admixes')

    def __call__(self, **kwargs):
        old_tree=kwargs['old_tree']
        return get_number_of_admixes(old_tree)
    
    def summary_of_phylogeny(self, tree):
        return get_number_of_admixes(tree)
    
    def make_histogram(self, data):
        return Summary.make_histogram(self, data, bins=range(10))
    
    def add_histogram(self, a):
        return Summary.add_histogram(self,  a, bins=range(10))

class s_branch_length(Summary):

    def __init__(self):
        super(s_branch_length,self).__init__('branch_length')

    def __call__(self, **kwargs):
        old_tree=kwargs['old_tree']
        return sum(get_all_branch_lengths(old_tree))
    
    def summary_of_phylogeny(self, tree):
        return sum(get_all_branch_lengths(tree))
    
class s_variable(Summary):
    
    def __init__(self, variable, pandable=True):
        super(s_variable, self).__init__(variable, pandable)

    def __call__(self, **kwargs):
        return kwargs[self.name]

 