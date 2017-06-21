from Rtree_operations import remove_parent_attachment, get_all_branch_descendants_and_rest, get_parent_of_branch, move_node, pretty_string, get_branch_length, update_node, get_all_branches, node_is_admixture
from nearest_branch_search import distanced_branch_lengths
from numpy.random import normal, choice
from copy import deepcopy

class sliding_rescale_class(object):
    
    new_nodes=1
    proposal_name='sliding_regraft'
    
    def __call__(self,*args, **kwargs):
        return make_sliding_rescale(*args, **kwargs)
    
def simulate_regraft_distance(param):
    return choice(3)+1

def get_thinned_pieces(tree,rescale_key, rescale_branch, cutoff_distance, parent_key):
    pieces= distanced_branch_lengths(tree, parent_key, visited_keys=[], upper_cap=cutoff_distance)
    thinned_pieces=[piece for piece in pieces if  piece.within_distance(cutoff_distance)]
    return thinned_pieces

def updater(param):
    def upd():
        return normal()*param

def resimulate_pieces(tree, pieces, max_distance, param, admixture_multiplier=1.0):
    visited_keys=[]
    upd=updater(param)
    for piece in pieces:
        if piece.child_key!='r':
            visited_keys.append(piece.parent_key)
            distance=piece.get_start_distance()
            d=(max_distance-distance)/max_distance
    
            if piece.child_branch==1:
                tree[piece.child_key][2]+= normal()*param*admixture_multiplier
            tree[piece.child_key][3+piece.child_branch]+= normal()*param
    return tree
                
        
    
def make_sliding_rescale(tree, param=0.1,pks={}):
    '''
    
    '''
    
    possible_branches= get_all_branches(tree)
        
    new_tree= deepcopy(tree)
    rescale_key, rescale_branch= possible_branches[choice(len(possible_branches), 1)[0]]
    pks['rescale_key']=rescale_key
    pks['rescale_branch']=rescale_branch
    pks['sliding_regraft_adap_param']= param
    
    cutoff_distance= simulate_regraft_distance(param)
    parent_key= get_parent_of_branch(tree, rescale_key, rescale_branch)
    thinned_pieces_forward=get_thinned_pieces(tree, rescale_key, rescale_branch, cutoff_distance, parent_key)
    
    print 'cutoff_distance', cutoff_distance
    print 'pieces:'
    for tp in thinned_pieces_forward:
        print tp
    
    new_tree=resimulate_pieces(new_tree, thinned_pieces_forward, cutoff_distance, param)

    
    forward=1.0
    backward=1.0
    
    return new_tree, forward, backward

if __name__=='__main__':
    from Rtree_operations import create_trivial_tree, pretty_string
    tree= create_trivial_tree(3,1.0)
    print pretty_string(tree)
    pks={}
    n,f,b= make_sliding_rescale(tree, pks=pks)
    print 
    print pretty_string(n)
    for e,v in pks.items():
        print e,':', v
    