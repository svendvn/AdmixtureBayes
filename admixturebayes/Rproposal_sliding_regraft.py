from Rproposal_regraft import get_possible_regrafters, thin_out_sibling
from Rtree_operations import remove_parent_attachment, get_all_branch_descendants_and_rest, get_parent_of_branch, move_node, pretty_string, get_branch_length, update_branch_length
from nearest_branch_search import distanced_branch_lengths
from scipy.stats import chi2
from numpy.random import choice
from copy import deepcopy
from random import getrandbits
from math import exp
import exponential_F
import gamma_restricted

class sliding_regraft_class(object):
    
    new_nodes=1
    proposal_name='sliding_regraft'
    adaption=True
    
    def __call__(self,*args, **kwargs):
        return make_sliding_regraft(*args, **kwargs)
    
class sliding_regraft_class_resimulate(object):
    new_nodes=1
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
    admixture_change=0
    proposal_name='sliding_regraft'
    adaption=True
    reverse='sliding_regraft'
    
    def __init__(self, param=False):
        self.param=param
    
    def __call__(self,*args, **kwargs):
        return make_sliding_regraft(*args, resimulate_moved_branch_length=self.param, **kwargs)
    
def simulate_regraft_distance(param):
    return chi2.rvs(1)*param

def get_thinned_pieces(tree,regraft_key, regraft_branch, distance_to_regraft, parent_key):
    pieces= distanced_branch_lengths(tree, parent_key, visited_keys=[regraft_key], upper_cap=distance_to_regraft)
    children, other= get_all_branch_descendants_and_rest(tree, regraft_key, regraft_branch)
    candidates= thin_out_sibling(tree, other, regraft_key)+[('r',0)]
    thinned_pieces=[piece for piece in pieces if (piece.get_branch_key() in candidates and piece.contains_distance(distance_to_regraft) )]
    return thinned_pieces
    
    
def make_sliding_regraft(tree, new_node=None, param=0.1, resimulate_moved_branch_length=False, pks={}):
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
    
    #print chosen_piece
    if new_node is None:
        new_node=str(getrandbits(68)).strip()
    new_tree =move_node(new_tree, regraft_key, regraft_branch, parent_key, distance_to_regraft, chosen_piece, new_node_name=new_node)
    
    parent_key= get_parent_of_branch(new_tree, regraft_key, regraft_branch)
    thinned_pieces_backward=get_thinned_pieces(new_tree, regraft_key, regraft_branch, distance_to_regraft, parent_key)
    
    backward_choices=len(thinned_pieces_backward)
    
    
    if resimulate_moved_branch_length:
        new_tree, forward, backward= resimulate_moved_branch(new_tree, regraft_key, regraft_branch, chosen_piece.get_lattitude(distance_to_regraft), resimulate_moved_branch_length)
    else:
        forward=1.0
        backward=1.0
    
    return new_tree, forward/forward_choices, backward/backward_choices

def resimulate_moved_branch(tree, key, branch, delta_L, alpha):
    old_length= get_branch_length(tree, key, branch)
    new_length= gamma_restricted.rvs(old_length, -delta_L, alpha)
    logpdf=gamma_restricted.logpdf(new_length, old_length, -delta_L, alpha)
    logpdf_back=gamma_restricted.logpdf(old_length, new_length, delta_L, alpha)
    backward=exp(-logpdf+logpdf_back) #for stability, both forward and backward are put in backward. 
    ratio= exp(logpdf-logpdf_back)
    update_branch_length(tree, key, branch, new_length)
    #print 'resimulated', old_length, 'to', new_length, 'on branch', str((key,branch)), 'and change in lattitude at', delta_L, 'backward',  backward
    return tree, 1.0, backward
    

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
    print pretty_string(t)
    letters=['n'+str(i) for i in xrange(100000)]
    for i in xrange(100):
        t,f,b= make_sliding_regraft(t, new_node=letters[i] , resimulate_moved_branch_length=True)
        print 'forward', f
        print 'backward',b
        print 'tree'
        print pretty_string(t)

    