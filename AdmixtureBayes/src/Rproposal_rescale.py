from Rtree_operations import update_all_branches, tree_on_the_border2_with_children
from copy import deepcopy
from numpy.random import normal

class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)

def rescale(tree, sigma=0.01, pks={}):
    new_tree=deepcopy(tree)
    updat=updater(sigma)
    return update_all_branches(new_tree, updat)


if __name__=='__main__':
    from tree_plotting import plot_graph
    plot_graph(tree_on_the_border2_with_children)
    new_tree=rescale(tree_on_the_border2_with_children)
    plot_graph(new_tree)
    
    print 'old_tree', tree_on_the_border2_with_children
    print 'new_tree', new_tree
    