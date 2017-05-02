from Rproposal_regraft import get_possible_regrafters, thin_out_sibling
from Rtree_operations import remove_parent_attachment, get_all_branch_descendants_and_rest, get_parent_of_branch, move_node, pretty_string
from nearest_branch_search import distanced_branch_lengths
from scipy.stats import chi2
from numpy.random import choice
from copy import deepcopy
from random import getrandbits

class sliding_regraft_class(object):
    
    new_nodes=1
    proposal_name='sliding_regraft'
    
    def __call__(self,*args, **kwargs):
        return make_sliding_regraft(*args, **kwargs)
    
def simulate_regraft_distance(param):
    return chi2.rvs(1)*param

def get_thinned_pieces(tree,regraft_key, regraft_branch, distance_to_regraft, parent_key):
    pieces= distanced_branch_lengths(tree, parent_key, visited_keys=[regraft_key], upper_cap=distance_to_regraft)
    children, other= get_all_branch_descendants_and_rest(tree, regraft_key, regraft_branch)
    candidates= thin_out_sibling(tree, other, regraft_key)+[('r',0)]
    thinned_pieces=[piece for piece in pieces if (piece.get_branch_key() in candidates and piece.contains_distance(distance_to_regraft) )]
    return thinned_pieces
    
    
def make_sliding_regraft(tree, new_node=None, param=0.1, pks={}):
    '''
    
    '''
    
    possible_nodes= get_possible_regrafters(tree)
        
    new_tree= deepcopy(tree)
    regraft_key, regraft_branch= possible_nodes[choice(len(possible_nodes), 1)[0]]
    pks['regraft_key']=regraft_key
    pks['regraft_branch']=regraft_branch
    pks['sliding_regraft_adap_param']= param
    
    distance_to_regraft= simulate_regraft_distance(param)
    parent_key= get_parent_of_branch(tree, regraft_key, regraft_branch)
    thinned_pieces_forward=get_thinned_pieces(tree, regraft_key, regraft_branch, distance_to_regraft, parent_key)
        
    forward_choices=len(thinned_pieces_forward)
    chosen_piece=thinned_pieces_forward[choice(len(thinned_pieces_forward),1)[0]]
    
    #print regraft_key, regraft_branch, chosen_piece
    if new_node is None:
        new_node=str(getrandbits(68)).strip()
    new_tree=move_node(new_tree, regraft_key, regraft_branch, parent_key, distance_to_regraft, chosen_piece, new_node_name=new_node)
    
    parent_key= get_parent_of_branch(new_tree, regraft_key, regraft_branch)
    thinned_pieces_backward=get_thinned_pieces(new_tree, regraft_key, regraft_branch, distance_to_regraft, parent_key)
    
    backward_choices=len(thinned_pieces_backward)
    
    return new_tree, 1.0/forward_choices, 1.0/backward_choices

if __name__=='__main__':
    
    test_tree={ 's3': ['x', None, None, 0.3, None, None, None],
                's2': ['a', None, None, 0.05, None, None, None],
                's1': ['d', None, None, 0.1, None, None, None],
                's4': ['b', None, None, 0.3, None, None, None],
                'a': ['b', 'c', 0.5, 0.2, 0.1, 's2', None],
                'c': ['x', 'd', 0.5, 0.083428092356333014, 0.1, 'a', None],
                'b': ['f', None, None, 0.05, None, 's4', 'a'],
                'd': ['r', None, None, 0.05, None, 's1', 'c'],
                'f': ['r', None, None, 0.02, None, 'b', 'x'],
                'x': ['f', None, None, 0.066571907643666994, None, 's3', 'c']}
    
    
    
    from Rcatalogue_of_trees import tree_good
    t=tree_good
    for _ in xrange(1000):
        t,f,b= make_sliding_regraft(t)
        print 'forward', f
        print 'backward',b
        print 'tree'
        print pretty_string(t)

    