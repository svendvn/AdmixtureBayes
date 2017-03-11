
tree_clean={'s1':['s1s2',None, None, 0.1,None],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3':['r',None, None, 0.2, None]}

tree_one_admixture={'s1':['s1b',None, None, 0.1,None],
      's1b':['s1s2','s3b',0.2, 0.1,0.2],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r',None, None, 0.2, None],
      's3':['s3b',None,None,0.2,None]}

tree_two_admixture={'s1':['s1b',None, None, 0.1,None],
      's1c':['s1s2','s3b', 0.4,0.05,0.1],
      's1b':['s1c','s3a',0.2, 0.05,0.2],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r',None, None, 0.2, None],
      's3':['s3a',None,None,0.1,None],
      's3a':['s3b', None,None,0.1,None]
      }

tree_two_admixture_cross={'s1':['s1b',None, None, 0.1,None],
      's1c':['s1s2','s3a', 0.4,0.05,0.1],
      's1b':['s1c',None,None, 0.05,None],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r','s1b', 0.4, 0.2, 0.2],
      's3':['s3a',None,None,0.1,None],
      's3a':['s3b', None,None,0.1,None]
      }

tree_illegal={'s1':['s1b',None, None, 0.1,None],
      's1c':['s1s2','s3a', 0.4,0.05,0.1],
      's1b':['s1c','s3b',0.2, 0.05,0.2],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r',None, None, 0.2, None],
      's3':['s3a',None,None,0.1,None],
      's3a':['s3b', None,None,0.1,None]
      }

tree_on_the_border={'s1':['c',None, None, 0.1,None],
      's2':['a',None, None,0.05,None],
      's3':['b',None,None, 0.3,None],
      'a':['b','d', 0.5,0.2,0.1],
      'c':['r','d',0.5,0.1,0.1],
      'd':['e',None,None,0.05,None],
      'b':['e',None,None,0.02,None],
      'e':['r',None,None,0.05,None]}

tree_on_the_border2={'s1':['d',None, None, 0.1,None],
      's2':['a',None, None,0.05,None],
      's3':['e',None,None, 0.3,None],
      's4':['b',None,None, 0.3,None],
      'a':['b','c', 0.5,0.2,0.1],
      'c':['e','d',0.5,0.1,0.1],
      'b':['f',None,None,0.05,None],
      'f':['r',None,None,0.02,None],
      'e':['f',None,None,0.05,None],
      'd':['r',None,None,0.05,None]}

tree_on_the_border2_with_children={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}


tree_admix_to_child={
    's1':['r',None,None, 0.1,None],
    's2':['s2a',None,None,0.1,None],
    's3':['s3s2',None,None,0.1,None],
    's2a':['s3s2','s3s2a', 0.5,0.1,0.13],
    's3s2':['s3s2a',None,None,0.1,None],
    's3s2a':['r',None,None,0.01]
    }

def create_trivial_tree(size, total_height=1.0):
    step_size=total_height/size
    tree={'s1':['n1',None,None,step_size,None, None,None],
          's2':['n1',None,None,step_size,None, None,None],
          'n1':['n2',None,None,step_size,None, 's1','s2']}
    nex_inner_node='n2'
    new_inner_node='n1'
    for k in range(3,size+1):
        old_inner_node='n'+str(k-2)
        new_inner_node='n'+str(k-1)
        nex_inner_node='n'+str(k)
        new_leaf='s'+str(k)
        tree[new_leaf]=[new_inner_node, None,None, step_size*(k-1),None, None,None]
        tree[new_inner_node]=[nex_inner_node, None,None, step_size, None, new_leaf,old_inner_node]
    del tree[new_inner_node]
    return _rename_root(tree, new_inner_node)

def extend_branch(node, pkey, grand_parent_key, p_to_gp):
    #print node, pkey, grand_parent_key, p_to_gp
    if node[0]==pkey:
        node[0]=grand_parent_key
        u=node[3]/(node[3]+p_to_gp)
        node[3]+=p_to_gp
        return node,u
    elif node[1]==pkey:
        node[1]=grand_parent_key
        u=node[4]/(node[4]+p_to_gp)
        node[4]+=p_to_gp
        return node,u
    else:
        assert False, 'extension of branch was not possible'

    
def remove_parent_attachment(tree, orphanota_key):
    '''
    This takes the tree and removes the parent of orphanonte_key.
    '''
    pkey=get_parents(tree[orphanota_key])[0]
    
    
    if pkey=='r':
        return remove_root_attachment(tree, orphanota_key)
    grand_pkey=get_parents(tree[pkey])[0]
    child_of_parent=get_other_children(tree[pkey], orphanota_key)[0]
    sib_node=tree[child_of_parent]
    tree[child_of_parent],u=extend_branch(sib_node, pkey, grand_pkey, tree[pkey][3])
    del tree[pkey]
    if grand_pkey!='r':
        tree[grand_pkey]=_rename_child(tree[grand_pkey], pkey, child_of_parent)
    tree[orphanota_key][0]=None
    return tree,"u",u

def remove_root_attachment(tree, orphanota_key):
    '''
    The situation is different when the root is removed because of the special naming strategy.
    Here a new root is born.
    '''
    rooted_keys=_find_rooted_nodes(tree)
    for key,len_to_root in rooted_keys:
        if key!=orphanota_key:
            if node_is_coalescence(tree[key]):
                tree=_rename_root(tree, key)
                r=len_to_root
                del tree[key]
            else:
                tree[key],r=get_branch_length_and_reset(tree[key], 'r', 'closed_branch')
                print 'closed_branch!'
            tree[orphanota_key][0]=None
    return tree,'r', r
    
def get_branch_length_and_reset(node, parent_key, new_length):
    if node[0]==parent_key:
        old_length=node[3]
        node[3]=new_length
        return node, old_length
    elif node[1]==parent_key:
        old_length=node[4]
        node[4]=new_length
        return node, old_length
    else:
        assert False, 'could not give new length because the parent was unknown'
    
def _rename_root(tree, old_name):
    for _,node in tree.items():
        if node[0]==old_name:
            node[0]='r'
        if (node[1] is not None and node[1]==old_name):
            node[1]='r'
    return tree

def _rename_child(node, old_name, new_name):
    if node[5]==old_name:
        node[5]=new_name
    elif node[6]==old_name:
        node[6]=new_name
    else:
        assert False, "tried to rename a child that didn't exist in its parents documents."
    return node
    
    
def _find_rooted_nodes(tree):
    res=[]
    for key,node in tree.items():
        if node[0]=='r' or (node[1] is not None and node[1]=='r'):
            if node[0]=='r':
                res.append((key,node[3]))
            else:
                res.append((key,node[4]))
    return res
    

def node_is_non_admixture(node):
    return (node[1] is None)

def node_is_admixture(node):
    return (node[1] is not None)

def node_is_coalescence(node):
    return (node[1] is None and node[5] is not None)

def node_is_leaf_node(node):
    return (node[1] is None and node[5] is None)



def get_descendants_and_rest(tree, key):
    all_keys=tree.keys()
    #print tree, key
    descendant_keys=_get_descendants(tree, key)
    return descendant_keys, list(set(all_keys)-set(descendant_keys))

def get_other_children(node, child_key):
    res=[]
    for n in get_children(node):
        if n is not None and n!=child_key:
            res.append(n)
    return res
    

def get_children(node):
    return node[5:7]

def _get_index_of_parent(node, parent):
    if node[0]==parent:
        return 0
    if node[1]==parent:
        return 1
    return -1

def has_child_admixture(tree, key):
    node=tree[key]
    child1,child2=get_children(node)
    if child1 is not None:
        if node_is_admixture(tree[child1]):
            return True
    if child2 is not None:
        return (node_is_admixture(tree[child2]))
    return False
        
def _get_descendants(tree, key):
    if tree[key][5] is None:
        return [key]
    else:
        ans=[key]+_get_descendants(tree, tree[key][5])
        if tree[key][6] is None:
            return ans
        return ans+_get_descendants(tree, tree[key][6])
    
def insert_children_in_tree(tree):
    children={key:[] for key in tree}
    for key in tree:
        parents = get_real_parents(tree[key])
        for parent in parents:
            if parent!='r':
                children[parent].append(key)
    #print children
    for key in tree:
        tree[key]=_update_parents(tree[key], children[key])
    return tree

def graft(tree, remove_key, add_to_branch, insertion_spot, new_node_code, which_branch):
    #if we are at the root, things are different.
    if add_to_branch=='r':
        return graft_onto_root(tree, insertion_spot, remove_key, new_name_for_old_root=new_node_code)
    
    #updating the grafted_branch. Easy.
    tree[remove_key][0]=new_node_code
    
    #dealing with the other child of the new node
    index_of_pop_name_to_change=which_branch
    index_of_branchl_to_change=which_branch+3
    #saving info:
    parent_of_branch=tree[add_to_branch][index_of_pop_name_to_change]
    length_of_piece_to_break_up=tree[add_to_branch][index_of_branchl_to_change]
    #updating node
    tree[add_to_branch][index_of_pop_name_to_change]=new_node_code
    tree[add_to_branch][index_of_branchl_to_change]=length_of_piece_to_break_up*insertion_spot
    
    #taking care of the new node
    tree[new_node_code]=[parent_of_branch, None, None, length_of_piece_to_break_up*(1.0-insertion_spot), None, add_to_branch, remove_key]

    #taking care of the parent of the new node
    if parent_of_branch != 'r':
        tree[parent_of_branch]=_update_child(tree[parent_of_branch], add_to_branch, new_node_code)
        
    return tree
    
def graft_onto_root(tree, insertion_spot, remove_key, new_name_for_old_root):    
    root_keys=_find_rooted_nodes(tree)

    if len(root_keys)==1:#this is the special case where an admixture leads to the new root
        return graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_keys[0])
    
    #updating the grafted_branch. Easy.
    tree[remove_key][0]='r'
    
    #dealing with the other child of the new node, but since the new node is the root, the old root is the new node. If that makes sense.
    #print 'root_keys', root_keys
    tree[new_name_for_old_root]=['r', None, None, insertion_spot, None,root_keys[0][0], root_keys[1][0]]
    
    #dealing with the children of the new node.
    for rkey, b_l in root_keys:
        tree[rkey]=_update_parent(tree[rkey], 'r', new_name_for_old_root)
    
    return tree

def graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_key):
    print 'undoing a closed branch', insertion_spot, remove_key, root_key
    tree[remove_key][0]='r'
    tree[root_key[0]],_=get_branch_length_and_reset(tree[root_key[0]], 'r', insertion_spot)
    return tree

def get_number_of_admixes(tree):
    return sum((1 for node in tree.values() if node_is_admixture(node)))
            


def _update_parents(node, new_parents):
    if len(new_parents)==1:
        res=node[:5]+[new_parents[0],None]
        return res
    if len(new_parents)==2:
        res=node[:5]+new_parents
        return res
    if len(new_parents)==0:
        res=node[:5]+[None]*2
        return res
    assert False, 'how many parents do you think you have?'
    
def _update_parent(node, old_parent, new_parent):
    if node[0]==old_parent:
        node[0]=new_parent
    elif node[1]==old_parent:
        node[1]=new_parent
    else:
        assert False, 'parent could not be updated'
    return node
        
def _update_child(node, old_child, new_child):
    if node[5]==old_child:
        node[5]=new_child
    elif node[6]==old_child:
        node[6]=new_child
    else:
        assert False, 'child could not be updated'
    return node
        
def get_parents(node):
    return node[:2]

def insert_admixture_node_halfly(tree, source):

def get_real_parents(node):
    ps=node[:2]
    return [p for p in ps if p is not None]

def get_other_parent(node, parent_key):
    if parent_key==node[0]:
        return node[1]
    elif parent_key==node[1]:
        return node[0]
    else:
        assert False, 'the shared parent was not a parent of the sibling.'
        
def halfbrother_is_uncle(tree, key):
    '''
    key is a non coalescence node. This function checks if a sibling is an admixture that goes to the parent and the grandparent at the same time.
    If removed, there would be a loop where one person has two of the same parent.
    '''
    parent_key=tree[key][0]
    sibling_key=get_other_children(tree[parent_key], key)[0]
    if node_is_non_admixture(tree[sibling_key]):
        return False
    bonus_parent=get_other_parent(tree[sibling_key], parent_key)
    grand_parent_key=tree[parent_key][0]
    return bonus_parent==grand_parent_key
        

def is_root(*keys):
    ad=[key=='r' for key in list(keys)]
    return any(ad)

def to_aarhus_admixture_graph(tree):
    leaves=[]
    inner_nodes=[]
    edges=[]
    admixture_proportions=[]
    for key,node in tree.items():
        if node_is_leaf_node(node):
            leaves.append([key])
        else:
            inner_nodes.append([key])
            if node_is_admixture(node):
                admixture_proportions.append([key, node[0],node[1],node[2]])
        ps=get_real_parents(node)
        for p in ps:
            edges.append([key,p])
    return leaves, inner_nodes, edges, admixture_proportions


    
if __name__=='__main__':
        tree={'s1':['s1s2',None, None, 0.1,None],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3':['r',None, None, 0.4, None]}
        
        from copy import deepcopy
        
        tree2=insert_children_in_tree(tree_on_the_border2)
        print tree2
        print tree_on_the_border2_with_children
        
        tree3=remove_parent_attachment(deepcopy(tree2), "s1")[0]
        print tree3
        print graft(tree3, 'f', 'r', 0.3, 'new_code', 0)
        #print remove_parent_attachment(tree2, 'f')
        #print get_descendants_and_rest(tree2, 'b')
        #print is_root(*get_parents(tree2['a']))
        
        trouble={'a': ['132', 'c', 0.5, 0.06013128348912011, 0.1, 's2', None], 'c': ['212', 'r', 0.5, 0.03639623632910125, 0.15000000000000002, 'a', None], 'b': ['r', None, None, 0.07, None, '132', 's4'], 'e': ['132', None, None, 0.05, None, '212', 's3'], 's3': ['e', None, None, 0.3, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['212', None, None, 0.1, None, None, None], 's4': ['b', None, None, 0.3, None, None, None], '132': ['b', None, None, 0.1398687165108799, None, 'a', 'e'], '212': ['e', None, None, 0.06360376367089875, None, 'c', 's1']}
        trouble= remove_parent_attachment(trouble, 'b')[0]
        print trouble
        print graft(trouble, 'b', 'r', 0.314, 'dont see me', 'hallo')
        
        trouble3={'a': ['n17', 'c', 0.5, 0.0006670327290825764, 0.1, 's2', None], 'c': ['n15', 'r', 0.5, 0.02087163982263861, 0.4814480657456043, 'a', None], 'n16': ['n17', None, None, 0.005272434567465561, None, 's4', 's3'], 'n17': [None, None, None, 0.013899593800954894, None, 'a', 'n16'], 'n15': ['r', None, None, 0.05969046586907494, None, 'c', 's1'], 's3': ['n16', None, None, 0.07815645814883887, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n15', None, None, 0.5947563021746359, None, None, None], 's4': ['n16', None, None, 0.00017898147838901196, None, None, None]}
        print graft(trouble3, 'n17', 'a', 1, 'n18', 1)
        
        tree_trouble={'a': ['n37', 'c', 0.5, 1.5717637721311875, 0.1, 's1', None], 'n66': ['r', None, None, 0.008798782728668674, None, 's3', 's4'], 'c': ['n54', 'n37', 0.5, 0.771318479326775, 0.07345113788460944, 'a', None], 's3': ['n66', None, None, 0.010969920361510089, None, None, None], 's2': ['n54', None, None, 0.404441491678861, None, None, None], 's1': ['a', None, None, 0.06451508173696463, None, None, None], 's4': ['n66', None, None, 1.7305330689019498, None, None, None], 'n67': ['r', None, None, 0.24519067463109384, None, 'n54', 'n37'], 'n54': ['n67', None, None, 0.25870104556004564, None, 'c', 's2'], 'n37': ['n67', None, None, 0.9342460567572629, None, 'c', 'a']}
        print 'tree_trouble', tree_trouble
        removed=remove_parent_attachment(tree_trouble, 'n54')[0]
        print 'pruned tree', removed
        adm=graft(removed, 'n54', 'c', 0.3, 'n68', 0) #FIXME: the function is being called like this via make_regraft
        print adm
        


        
        
        