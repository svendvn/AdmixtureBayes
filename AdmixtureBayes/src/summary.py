import seaborn as sns

def base_trajectory(facet_grid, summary_name):
    return facet_grid.map(sns.tsplot, summary_name)


def base_histogram(facet_grid, summary_name, **kwargs):
    return facet_grid.map(sns.distplot, summary_name, kde=False, **kwargs)

#from Rtree_operations import get_number_of_admixes, get_all_branch_lengths

class Summary(object):
       
    def __init__(self, name, pandable=True):
        self.name=name
        self.pandable=pandable
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return str(output)
    
    def make_trajectory(self, facet_grid):
        return base_trajectory(facet_grid, self.name)
    
    def make_histogram(self, facet_grid, **kwargs):
        return base_histogram(facet_grid, self.name, **kwargs)
    
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

 