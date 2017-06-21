from Rtree_operations import remove_parent_attachment, get_all_branch_descendants_and_rest, get_parent_of_branch, move_node, pretty_string, get_branch_length, update_node, get_all_branches, node_is_admixture
from nearest_branch_search import distanced_branch_lengths
from numpy.random import normal
from boto.roboto import param

class sliding_regraft_class(object):
    
    new_nodes=1
    proposal_name='sliding_regraft'
    
    def __call__(self,*args, **kwargs):
        return make_sliding_regraft(*args, **kwargs)
    
class sliding_regraft_class_resimulate(object):
    new_nodes=1
    proposal_name='sliding_regraft'
    
    def __init__(self, param):
        self.param=param
    
    def __call__(self,*args, **kwargs):
        return make_sliding_regraft(*args, resimulate_moved_branch_length=self.param, **kwargs)
    
def simulate_regraft_distance(param):
    return chi2.rvs(1)*param

def get_thinned_pieces(tree,regraft_key, regraft_branch, cutoff_distance, parent_key):
    pieces= distanced_branch_lengths(tree, parent_key, visited_keys=[], upper_cap=cutoff_distance)
    children, other= get_all_branch_descendants_and_rest(tree, regraft_key, regraft_branch)
    candidates= thin_out_sibling(tree, other, regraft_key)+[('r',0)]
    thinned_pieces=[piece for piece in pieces if  piece.within_distance(cutoff_distance)]
    return thinned_pieces

def updater(param):
    def upd():
        return normal()*param

def resimulate_pieces(tree, pieces, max_distance, param, admixture_multiplier=1.0):
    visited_keys=[]
    upd=updater(param)
    for piece in pieces:
        visited_keys.append(piece.parent_key)
        distance=piece.get_distance()
        d=(max_distance-distance)/max_distance

        if piece.child_branch==1:
            tree[piece.child_key][2]+= normal()*param*admixture_multiplier
        tree[piece.child_key][3+piece.child_key]+= normal()*param
    return param
                
        
    
def make_sliding_regraft(tree, new_node=None, param=0.1, resimulate_moved_branch_length=False, pks={}):
    '''
    
    '''
    
    possible_branches= get_all_branches(tree)
        
    new_tree= deepcopy(tree)
    rescale_key, rescale_branch= possible_nodes[choice(len(possible_branches), 1)[0]]
    pks['regraft_key']=regraft_key
    pks['regraft_branch']=regraft_branch
    pks['sliding_regraft_adap_param']= param
    
    cutoff_distance= simulate_regraft_distance(param)
    parent_key= get_parent_of_branch(tree, regraft_key, regraft_branch)
    thinned_pieces_forward=get_thinned_pieces(tree, regraft_key, regraft_branch, cutoff_distance, parent_key)
        
    new_tree=resimulate_pieces(new_tree, thinned_pieces_forward, cutoff_distance, param)
    
    chosen_piece=thinned_pieces_forward[choice(len(thinned_pieces_forward),1)[0]]
    
    #print chosen_piece
    if new_node is None:
        new_node=str(getrandbits(68)).strip()
    new_tree =move_node(new_tree, regraft_key, regraft_branch, parent_key, cutoff_distance, chosen_piece, new_node_name=new_node)
    
    parent_key= get_parent_of_branch(new_tree, regraft_key, regraft_branch)
    thinned_pieces_backward=get_thinned_pieces(new_tree, regraft_key, regraft_branch, cutoff_distance, parent_key)
    
    backward_choices=len(thinned_pieces_backward)
    
    
    if resimulate_moved_branch_length:
        new_tree, forward, backward= resimulate_moved_branch(new_tree, regraft_key, regraft_branch, chosen_piece.get_lattitude(cutoff_distance), resimulate_moved_branch_length)
    else:
        forward=1.0
        backward=1.0
    
    return new_tree, forward/forward_choices, backward/backward_choices

if __name__=='__main__':
    