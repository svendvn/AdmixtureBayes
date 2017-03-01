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
    u=sib_node[3]+tree[pkey][3]
    tree[child_of_parent]=[grand_pkey, None,None, u, None, sib_node[5],sib_node[6]]
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
            tree=_rename_root(tree, key)
            r=len_to_root
            del tree[key]
            tree[orphanota_key][0]=None
    return tree,'r', r
    
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
    print children
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
    
    #updating the grafted_branch. Easy.
    tree[remove_key][0]='r'
    
    #dealing with the other child of the new node, but since the new node is the root, the old root is the new node. If that makes sense.
    tree[new_name_for_old_root]=['r', None, None, insertion_spot, None]+root_keys
    
    #dealing with the children of the new node.
    for rkey in root_keys:
        tree[rkey]=_update_parent(tree[rkey], 'r', new_name_for_old_root)
    
    return tree


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

def get_real_parents(node):
    ps=node[:2]
    return [p for p in ps if p is not None]

def is_root(*keys):
    ad=[key=='r' for key in list(keys)]
    return any(ad)
    
if __name__=='__main__':
        tree={'s1':['s1s2',None, None, 0.1,None],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3':['r',None, None, 0.4, None]}
        
        
        tree2=insert_children_in_tree(tree_on_the_border2)
        print tree2
        print tree_on_the_border2_with_children
        
        print remove_parent_attachment(tree2, "b")
        #print get_descendants_and_rest(tree2, 'b')
        #print is_root(*get_parents(tree2['a']))
        