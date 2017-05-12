def create_trivial_tree(size, total_height=1.0):
    '''
    constructs tree of the form (..((s1,s2),s3),s4)...)
    '''
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
    return rename_root(tree, new_inner_node)

def create_trivial_equibranched_tree(size, height=1.0):
    '''
    constructs tree of the form (..((s1,s2),s3),s4)...)
    '''
    tree={'s1':['n1',None,None,height,None, None,None],
          's2':['n1',None,None,height,None, None,None],
          'n1':['n2',None,None,height,None, 's1','s2']}
    nex_inner_node='n2'
    new_inner_node='n1'
    for k in range(3,size+1):
        old_inner_node='n'+str(k-2)
        new_inner_node='n'+str(k-1)
        nex_inner_node='n'+str(k  )
        new_leaf='s'+str(k)
        tree[new_leaf]=[new_inner_node, None,None, height,None, None,None]
        tree[new_inner_node]=[nex_inner_node, None,None, height, None, new_leaf,old_inner_node]
    del tree[new_inner_node]
    return rename_root(tree, new_inner_node)

def find_children(tree, parent_key):
    res=[]
    for key,node in tree.items():
        if node[0]==parent_key:
            res.append(key)
        if node[1]==parent_key:
            res.append(key)
    while len(res)<2:
        res.append(None)
    return res

def create_balanced_tree(size, height=1.0):
    return finish_tree_with_coalescences({}, get_trivial_nodes(size), height)
        

def finish_tree_with_coalescences(tree, keys_to_finish, height=1.0):
    while len(keys_to_finish)>1:
        key1,key2=keys_to_finish[:2]
        if len(keys_to_finish)==2:
            new_key='r'
        else:
            new_key=key1+key2
        tree[key1]=[new_key,None,None,height,None]+find_children(tree,key1)
        tree[key2]=[new_key,None,None,height,None]+find_children(tree,key2)
        keys_to_finish=keys_to_finish[2:]
        keys_to_finish.append(new_key)
        print len(keys_to_finish)
    return tree
        

def create_burled_leaved_tree(size, height):
    res={}
    for i in range(size):
        res['s'+str(i+1)]=['a'+str(i+1), None, None, height, None,None,None]
        res['a'+str(i+1)]=['b'+str(i+1), 'm'+str(i+1), 0.5, height, height,'s'+str(i+1),None]
        res['b'+str(i+1)]=['n'+str(i+1), 'm'+str(i+1), 0.5, height, height,'a'+str(i+1),None]
        res['m'+str(i+1)]=['n'+str(i+1),None ,None, height, None,'a'+str(i+1),'b'+str(i+1)]
    return finish_tree_with_coalescences(res, ['n'+str(i+1) for i in range(size)], height)
    

def get_trivial_nodes(size):
    return ['s'+str(n+1) for n in xrange(size)]

def get_distance_to_root(tree, key, function=max):
    if key=='r':
        return 0.0
    node=tree[key]
    if node_is_admixture(node):
        return function(get_distance_to_root(tree, node[0])+node[3], 
                   get_distance_to_root(tree, node[1])+node[4], node[2])
    else:
        return get_distance_to_root(tree, node[0])+node[3]

def special_max(x,y,z=None):
    return max(x,y)

def special_min(x,y,z=None):
    return min(x,y)

def average_admixture_node(x,y,z):
    return z*x+(1-z)*y

def get_max_distance_to_root(tree):
    return max(get_leaf_distances_to_root(tree, function=special_max))

def get_min_distance_to_root(tree):
    return min(get_leaf_distances_to_root(tree, function=special_min))

def get_average_distance_to_root(tree):
    av_lengths=get_leaf_distances_to_root(tree, function=average_admixture_node)
    return float(sum(av_lengths))/len(av_lengths)

def get_leaf_distances_to_root(tree, function=max):
    res=[]
    for key, node in tree.items():
        if node_is_leaf_node(node):
            res.append(get_distance_to_root(tree, key, function=function))
    return res

def get_number_of_ghost_populations(tree):
    count=0
    for key,node in tree.items():
        if node_is_coalescence(node):
            child1,child2=get_children(node)
            if node_is_admixture(tree[child1]) and node_is_admixture(tree[child2]):
                count+=1
    return count

    
def update_all_branches(tree, updater):
    for key, node in tree.items():
        if node_is_admixture(node):
            node[2]+=updater()
            node[3]+=updater()
            node[4]+=updater()
            if node[2]<0 or node[2]>1 or node[3]<0 or node[4]<0:
                return None
        else:
            node[3]+=updater()
            if node[3]<0:
                return None        
    return tree

def update_branch_length(tree,key,branch, new_length):
    tree[key][branch+3]=new_length

def extend_branch(node, pkey, grand_parent_key, p_to_gp):
    #print node, pkey, grand_parent_key, p_to_gp
    if node[0]==pkey:
        node[0]=grand_parent_key
        u=node[3]/(node[3]+p_to_gp)
        node[3]+=p_to_gp
        return node,u,node[3]
    elif node[1]==pkey:
        node[1]=grand_parent_key
        u=node[4]/(node[4]+p_to_gp)
        node[4]+=p_to_gp
        return node,u,node[4]
    else:
        assert False, 'extension of branch was not possible'

    
def remove_parent_attachment(tree, orphanota_key, orphanota_branch):
    '''
    This takes the tree and removes the parent of orphanonte_key.
    '''
    pkey=get_parent_of_branch(tree, orphanota_key, orphanota_branch)
    
    if pkey=='r':
        return remove_root_attachment(tree, orphanota_key, orphanota_branch)
    grand_pkey=get_parents(tree[pkey])[0]
    child_of_parent=get_other_children(tree[pkey], orphanota_key)[0]
    sib_node=tree[child_of_parent]
    tree[child_of_parent],u,extended_branch_length=extend_branch(sib_node, pkey, grand_pkey, tree[pkey][3])
    del tree[pkey]
    if grand_pkey!='r':
        tree[grand_pkey]=_rename_child(tree[grand_pkey], pkey, child_of_parent)
    tree[orphanota_key][orphanota_branch]=None
    return tree,"u",u*extended_branch_length,extended_branch_length

def remove_root_attachment(tree, orphanota_key, orphanota_branch):
    '''
    The situation is different when the root is removed because of the special naming strategy.
    
                r
               / \
             /    \
      orphanota   new_root 
    Here a new root is born.
    
    '''
    rooted_keys=find_rooted_nodes(tree)
    for key,branch,len_to_root in rooted_keys:
        if key!=orphanota_key or orphanota_branch!=branch:
            if node_is_coalescence(tree[key]):
                tree=rename_root(tree, key)
                r=len_to_root
                del tree[key]
            else:
                tree[key],r=get_branch_length_and_reset(tree[key], 'r', 'closed_branch')
                #print 'closed_branch!'
            tree[orphanota_key][orphanota_branch]=None
    return tree,'r', r,None
    
def get_branch_length_and_reset(node, parent_key, new_length,add=False):
    if node[0]==parent_key:
        old_length=node[3]
        if add:
            node[3]+=new_length
        else:
            node[3]=new_length
        return node, old_length
    elif node[1]==parent_key:
        old_length=node[4]
        if add:
            node[4]+=new_length
        else:
            node[4]=new_length
        return node, old_length
    else:
        assert False, 'could not give new length because the parent was unknown'
    
def rename_root(tree, old_name):
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

def rename_parent(node, old_name, new_name):
    if node[0]==old_name:
        node[0]=new_name
    elif node[1]==old_name:
        node[1]=new_name
    else:
        assert False, 'tried to rename a parent that didnt exist in its childs document'
    return node
    
def find_rooted_nodes(tree):
    res=[]
    for key,node in tree.items():
        if node[0]=='r' or (node[1] is not None and node[1]=='r'):
            if node[0]=='r':
                res.append((key,0,node[3]))
            else:
                res.append((key,1,node[4]))
    return res

def _find_rooted_branches(tree):
    res=[]
    for key,node in tree.items():
        if node[0]=='r':
            res.append((key,0))
        if node[1] is not None and node[1]=='r':
            res.append((key,1))
    return res

def _get_root_sibling(tree, child_key, child_branch):
    root_keys=_find_rooted_branches(tree)
    if len(root_keys)==1:
        return child_key, child_branch
    else:
        if root_keys[0][0]==child_key and root_keys[0][1]==child_branch:
            return root_keys[1]
        else:
            return root_keys[0]

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

def _get_descendants(tree, key):
    if tree[key][5] is None:
        return [key]
    else:
        ans=[key]+_get_descendants(tree, tree[key][5])
        if tree[key][6] is None:
            return ans
        return ans+_get_descendants(tree, tree[key][6])

def get_all_branches(tree):
    res=[]
    for key, node in tree.items():
        if node_is_admixture(node):
            res.extend([(key, 0),(key,1)])
        else:
            res.append((key,0))
    return res

def get_all_branch_descendants_and_rest(tree, key,branch):
    all_branches=get_all_branches(tree)
    descendant_branches=_get_descendant_branches(tree, key,branch)
    return descendant_branches, list(set(all_branches)-set(descendant_branches))

def _get_descendant_branches(tree, key, branch):
    child_key1=tree[key][5]
    if child_key1 is None: #branch is leaf
        return [(key,0)]
    else:
        child_key1_branch=mother_or_father(tree, child_key1, key)
        ans=[(key, branch)]+_get_descendant_branches(tree, child_key1, child_key1_branch)
        child_key2=tree[key][6]
        if child_key2 is None:
            return ans
        return ans+_get_descendant_branches(tree, child_key2, mother_or_father(tree,child_key2, key))

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
        

    
def get_leaf_keys(tree):
    res=[]
    for key, node in tree.items():
        if node_is_leaf_node(node):
            res.append(key)
    return res

def get_no_leaves(tree):
    return len(get_leaf_keys(tree))

def get_categories(tree):
    leaves=[]
    admixture_nodes=[]
    coalescence_nodes=[]
    for key, node in tree.items():
        if node_is_leaf_node(node):
            leaves.append(key)
        if node_is_coalescence(node):
            coalescence_nodes.append(key)
        if node_is_admixture(node):
            admixture_nodes.append(key)
    return leaves, coalescence_nodes, admixture_nodes

def get_parent_of_branch(tree, key, branch):
    assert key!='r', 'Tried to access the parent of the root branch'
    return tree[key][branch]

def get_branch_length(tree,key,branch):
    assert key!='r', 'Tried to access the length of the root branch'
    return tree[key][branch+3]

def get_branch_length_from_parent(tree, child_key, parent_key):
    return tree[child_key][3+mother_or_father(tree, child_key, parent_key)]

def get_admixture_proportion(tree, child_key,child_branch):
    key=get_parent_of_branch(tree, child_key,child_branch)
    assert node_is_admixture(tree[key]), 'Tried to get the admixture proportion of a non-admixture node'
    return tree[key][2]

def update_parent_and_branch_length(tree, child_key, child_branch, new_parent, new_branch_length):
    assert child_key!='r', 'Tried to update the root branch'
    tree[child_key][child_branch]=new_parent
    tree[child_key][child_branch+3]=new_branch_length
    return tree
    


def get_destination_of_lineages(tree, ready_lineages):
    single_coalescences={} #list of tuples ((key,branch),(sister_key,sister_branch))
    double_coalescences=[]
    admixtures=[]
    for key, branch in ready_lineages:
        if (key,branch) in single_coalescences:
            double_coalescences.append(((key,branch),single_coalescences[(key,branch)]))
            del single_coalescences[(key,branch)]
            continue
        parent_key=tree[key][branch]
        if parent_key=='r':
            sister_key, sister_branch=_get_root_sibling(tree, key, branch)
            single_coalescences[(sister_key,sister_branch)]=(key,branch)
            continue
        parent=tree[parent_key]
        if node_is_coalescence(parent):
            sister_key, sister_branch=get_sister_branch(tree, parent,key, branch)
            single_coalescences[(sister_key,sister_branch)]=(key,branch)
        elif node_is_admixture(parent):
            admixtures.append((key,branch))
        else:
            assert False, 'the parent of a node was neither admixture nor coalescence'
    return double_coalescences, single_coalescences, admixtures

def pretty_string(tree):
    keys,vals=tree.keys(),tree.values()
    res=''
    res+='{ '+'\n'
    for key,val in zip(keys,vals):
        if node_is_leaf_node(val):
            res+='  '+key+': '+str(val)+'\n'
    res+='  ,'+'\n'
    for key,val in zip(keys,vals):
        if node_is_admixture(val):
            res+='  '+key+': '+str(val)+'\n'
    res+='  ,'+'\n'
    for key,val in zip(keys,vals):
        if node_is_coalescence(val):
            res+='  '+key+': '+str(val)+'\n'
    res+='}'
    return res

def pretty_print(tree):
    print pretty_string(tree)

def get_sister_branch(tree, parent, key, branch):
    #print parent,key
    #print pretty_string(tree)
    if parent[5]==parent[6]:
        return key, other_branch(branch)
    else:
        if parent[5]==key:
            return parent[6], mother_or_father(tree, parent[6], tree[key][branch])
        elif parent[6]==key:
            return parent[5], mother_or_father(tree, parent[5], tree[key][branch])
        else:
            assert False, "the parent was not really a parent"+'\n'+pretty_string(tree)+'\n'+'parent,key,branch='+str(parent)+','+str(key)+','+str(branch)
        
        

def propagate_married(tree, list_of_pairs):
    res=[]
    for (key1,branch1),(_,_) in list_of_pairs:
        parent_key=get_parent_of_branch(tree, key1, branch1)
        res.append((parent_key,0))
    return res

def propagate_admixtures(tree, list_of_admixtures):
    res=[]
    for key,branch in list_of_admixtures:
        parent_key=get_parent_of_branch(tree, key, branch)
        res.append((parent_key,0))
        res.append((parent_key,1))
    return res
    
        
def mother_or_father(tree, child_key, parent_key):
    if tree[child_key][0]==parent_key:
        return 0
    elif tree[child_key][1]==parent_key:
        return 1
    assert False, 'The child did not match its parent'+\
                  '\n'+pretty_string(tree)+'\n\n'+\
                  'child_key,parent_key='+str(child_key)+\
                  ','+str(parent_key)
    
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



def get_all_branch_lengths(tree):
    res=[]
    for key, node in tree.items():
        if node_is_non_admixture(node):
            res.append(node[3])
        else:
            res.extend(node[3:5])
    return res
    
def get_all_admixture_proportions(tree):
    res=[]
    for key, node in tree.items():
        if node_is_admixture(node):
            res.append(node[2])
    return res

def convert_to_vector(tree, keys=None):
    if keys is None:
        keys=tree.keys()
    res=[]
    for key in keys:
        node=tree[key]
        if node_is_non_admixture(node):
            res.append(node[3])
        else:
            res.extend(node[3:5])
            res.append(node[2])
    return res

def graft(tree, remove_key, add_to_branch, insertion_spot, new_node_code, which_branch, remove_branch=0):
    #if we are at the root, things are different.
    if add_to_branch=='r':
        return graft_onto_root(tree, insertion_spot, remove_key, new_name_for_old_root=new_node_code, remove_branch=remove_branch)
    
    #updating the grafted_branch. Easy.
    tree[remove_key][remove_branch]=new_node_code
    
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
    
def graft_onto_root(tree, insertion_spot, remove_key, new_name_for_old_root, remove_branch=0):    
    root_keys=find_rooted_nodes(tree)

    if len(root_keys)==1:#this is the special case where an admixture leads to the new root
        return graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_keys[0], remove_branch=remove_branch)
    
    #updating the grafted_branch. Easy.
    tree[remove_key][remove_branch]='r'
    
    #dealing with the other child of the new node, but since the new node is the root, the old root is the new node. If that makes sense.
    #print 'root_keys', root_keys
    tree[new_name_for_old_root]=['r', None, None, insertion_spot, None,root_keys[0][0], root_keys[1][0]]
    
    #dealing with the children of the new node.
    for rkey, _, b_l in root_keys:
        tree[rkey]=_update_parent(tree[rkey], 'r', new_name_for_old_root)
    
    return tree

def materialize_root(tree, new_key):
    (child_key1, child_branch1,_),(child_key2, child_branch2,_) = find_rooted_nodes(tree)
    tree[new_key]=['r',None,None,None,None,child_key1, child_key2]
    tree[child_key1][child_branch1]=new_key
    tree[child_key2][child_branch2]=new_key
    return tree

def move_node(tree, regraft_key, regraft_branch, parent_key, distance_to_regraft, chosen_piece, new_node_name='x'):
    if parent_key=='r':
        sister_key,sister_branch= _get_root_sibling(tree, regraft_key, regraft_branch)
        if chosen_piece.child_key=='r':
            u1,_=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            tree[sister_key][sister_branch+3]=tree[sister_key][sister_branch+3]+u1
        elif chosen_piece.child_key==sister_key and chosen_piece.child_branch==sister_branch:
            u1,_=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            tree[sister_key][sister_branch+3]=u1
        else:
            tree=rename_root(tree, sister_key)
            del tree[sister_key]
            u1,u2 =chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            if chosen_piece.parent_key == sister_key:
                tree[new_node_name]=['r', None,None,u2,None,regraft_key,chosen_piece.child_key]
            else:
                tree[new_node_name]=[chosen_piece.parent_key, None,None,u2,None,regraft_key,chosen_piece.child_key]
            if chosen_piece.parent_key=='r':
                pass
            elif chosen_piece.parent_key != sister_key:
                tree[chosen_piece.parent_key]=_rename_child(tree[chosen_piece.parent_key], chosen_piece.child_key, new_node_name)
            tree=update_parent_and_branch_length(tree, chosen_piece.child_key, chosen_piece.child_branch, new_node_name, u1)
            tree[regraft_key][regraft_branch]=new_node_name
    else:
        if chosen_piece.child_key=='r':
            tree=materialize_root(tree, new_node_name)
            u1,_=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            tree[new_node_name][3]=u1       
            tree[regraft_key][regraft_branch]='r'     
        else:
            u1,u2=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            #print u1,u2
            tree[new_node_name]=[chosen_piece.parent_key, None,None,u2,None,regraft_key,chosen_piece.child_key]
            if chosen_piece.parent_key=='r':
                pass
            else:
                tree[chosen_piece.parent_key]=_rename_child(tree[chosen_piece.parent_key], chosen_piece.child_key, new_node_name)
            tree=update_parent_and_branch_length(tree, chosen_piece.child_key, chosen_piece.child_branch, new_node_name, u1)
            tree[regraft_key][regraft_branch]=new_node_name
        sibling_key=get_other_children(tree[parent_key], regraft_key)[0]
        grandfather_key=tree[parent_key][0]
        _=get_branch_length_and_reset(tree[sibling_key],parent_key, tree[parent_key][3], add=True)
        tree[sibling_key]=rename_parent(tree[sibling_key], parent_key, grandfather_key)
        if grandfather_key!='r':
            tree[grandfather_key]=_rename_child(tree[grandfather_key], parent_key, sibling_key)
        del tree[parent_key]
    
    return tree
        
    

def graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_key, remove_branch=0):
    #print 'undoing a closed branch', insertion_spot, remove_key, root_key
    tree[remove_key][remove_branch]='r'
    tree[root_key[0]],_=get_branch_length_and_reset(tree[root_key[0]], 'r', insertion_spot)
    return tree

def change_admixture(node):
    assert node_is_admixture(node), 'tried to change admixture branches for a node without admixture'
    new_node=[node[1],node[0],1.0-node[2],node[4],node[3]]+node[5:]
    return new_node

def get_number_of_admixes(tree):
    return sum((1 for node in tree.values() if node_is_admixture(node)))

def get_number_of_leaves(tree):
    return sum((1 for node in tree.values() if node_is_leaf_node(node)))
            
def remove_admix(tree, rkey, rbranch):
    '''
    removes an admixture. besides the smaller tree, also the disappeared branch lengths are returned.
    parent_key          sparent_key
        |                |
        |t_1             | t_4
        |   __---- source_key
      rkey/   t_5       \
        |                \t_3
        |t_2          sorphanota_key  
    orphanota_key   
    and alpha=admixture proportion. The function returns new_tree, (t1,t2,t3,t4,t5), alpha 

    Note that t_4 could be None if source key is root. The source_key node is not an admixture node by assumption.
    '''
    rnode= tree[rkey]
    orphanota_key= get_children(rnode)[0]
    parent_key= rnode[other_branch(rbranch)]
    source_key= rnode[rbranch]
    t1= rnode[3+other_branch(rbranch)]
    t5= rnode[3+rbranch]
    alpha=rnode[2]
    
    tree[orphanota_key],t2=get_branch_length_and_reset(tree[orphanota_key], rkey, t1, add=True)
    tree[orphanota_key]=_update_parent(tree[orphanota_key], rkey, parent_key)
    orphanota_branch=_get_index_of_parent(tree[orphanota_key], parent_key)
    
    if parent_key!='r':
        tree[parent_key]=_update_child(tree[parent_key], old_child=rkey, new_child=orphanota_key)
    if source_key=='r':
        tree,_,t3,_=remove_root_attachment(tree, rkey) #now sorphanota_key is new root
        del tree[rkey]
        return tree, (t1,t2,t3,None,t5),alpha,(orphanota_key,orphanota_branch)
    del tree[rkey]
    source_node=tree[source_key]
    sorphanota_key=get_other_children(source_node, child_key=rkey)[0]
    sparent_key=source_node[0]
    t4=source_node[3]
    if sparent_key!='r':
        tree[sparent_key]=_update_child(tree[sparent_key], source_key, sorphanota_key)
    tree[sorphanota_key],t3=get_branch_length_and_reset(tree[sorphanota_key], source_key, t4, add=True)
    tree[sorphanota_key]=_update_parent(tree[sorphanota_key], source_key, sparent_key)
    del tree[source_key]
    return tree, (t1,t2,t3,t4,t5), alpha, (orphanota_key,orphanota_branch)

def remove_admix2(tree, rkey, rbranch, pks={}):
    '''
    removes an admixture. besides the smaller tree, also the disappeared branch lengths are returned.
    parent_key          sparent_key
        |                |
        |t_1             | t_4
        |   __---- source_key
      rkey/   t_5       \
        |                \t_3
        |t_2          sorphanota_key  
    orphanota_key   
    and alpha=admixture proportion. The function returns new_tree, (t1,t2,t3,t4,t5), alpha 

    Note that t_4 could be None if source key is root. The source_key node is not an admixture node by assumption.
    '''
    rnode= tree[rkey]
    orphanota_key= get_children(rnode)[0]
    parent_key= rnode[other_branch(rbranch)]
    source_key= rnode[rbranch]
    t1= rnode[3+other_branch(rbranch)]
    t5= rnode[3+rbranch]
    alpha=rnode[2]
    
    tree[orphanota_key],t2=get_branch_length_and_reset(tree[orphanota_key], rkey, t1, add=True)
    tree[orphanota_key]=_update_parent(tree[orphanota_key], rkey, parent_key)
    orphanota_branch=_get_index_of_parent(tree[orphanota_key], parent_key)
    pks['orphanota_key']=orphanota_key
    pks['orphanota_branch']=orphanota_branch
    
    if parent_key!='r':
        tree[parent_key]=_update_child(tree[parent_key], old_child=rkey, new_child=orphanota_key)
    if source_key=='r':
        tree,_,t3,_=remove_root_attachment(tree, rkey, 0) #now sorphanota_key is new root
        del tree[rkey]
        pks['sorphanota_key']='r'
        pks['sorphanota_branch']=0
        return tree, (t1,t2,t3,None,t5),alpha
    del tree[rkey]
    source_node=tree[source_key]
    sorphanota_key=get_other_children(source_node, child_key=rkey)[0]
    sparent_key=source_node[0]
    t4=source_node[3]
    if sparent_key!='r':
        tree[sparent_key]=_update_child(tree[sparent_key], source_key, sorphanota_key)
    tree[sorphanota_key],t3=get_branch_length_and_reset(tree[sorphanota_key], source_key, t4, add=True)
    tree[sorphanota_key]=_update_parent(tree[sorphanota_key], source_key, sparent_key)
    pks['sorphanota_key']=sorphanota_key
    pks['sorphanota_branch']=mother_or_father(tree, sorphanota_key, sparent_key)
    del tree[source_key]
    return tree, (t1,t2,t3,t4,t5), alpha
    
def other_branch(branch):
    if branch==0:
        return 1
    elif branch==1:
        return 0
    else:
        assert False, 'illegal branch'
    

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

def insert_admixture_node_halfly(tree, source_key, source_branch, insertion_spot, admix_b_length, new_node_name, admixture_proportion=0.51):
    node=tree[source_key]
    old_parent=node[source_branch]
    old_branch_length=node[source_branch+3]
    node[source_branch]=new_node_name
    node[source_branch+3]=old_branch_length*insertion_spot
    tree[new_node_name]=[old_parent, None, admixture_proportion, old_branch_length*(1-insertion_spot), admix_b_length, source_key, None]
    if old_parent != 'r':
        tree[old_parent]=_update_child(tree[old_parent], old_child=source_key, new_child=new_node_name)
    return tree

def get_real_parents(node):
    ps=node[:2]
    return [p for p in ps if p is not None]

def get_real_children(node):
    cs=node[5:7]
    return [c for c in cs if c is not None]

def get_real_children_root(tree, key):
    if key=='r':
        a,b=find_rooted_nodes(tree)
        return [(a[0],a[1]), (b[0],b[1])]
    else:
        c_keys=get_real_children(tree[key])
        res=[]
        for c_key in c_keys:
            res.append((c_key, mother_or_father(tree, c_key, key)))
        return res
        

def get_other_parent(node, parent_key):
    if parent_key==node[0]:
        return node[1]
    elif parent_key==node[1]:
        return node[0]
    else:
        assert False, 'the shared parent was not a parent of the sibling.'
        
def halfbrother_is_uncle(tree, key, parent_key):
    '''
    parent_key is a non admixture node. This function checks if a sibling is an admixture that goes to the parent and the grandparent at the same time.
    If removed, there would be a loop where one person has two of the same parent.
    '''
    sibling_key=get_other_children(tree[parent_key], key)[0]
    if node_is_non_admixture(tree[sibling_key]):
        return False
    bonus_parent=get_other_parent(tree[sibling_key], parent_key)
    grand_parent_key=tree[parent_key][0]
    return bonus_parent==grand_parent_key

def parent_is_spouse(tree, key, direction):
    '''
    key is an admixture node. This function checks if the parent in the direction of 'direction' also has a child with the admixture node.
    '''
    parent_key=tree[key][direction]
    offspring_key=get_real_children(tree[key])[0] #there is only one because key is an admixture node
    if node_is_non_admixture(tree[offspring_key]):
        return False
    spouse_key=get_other_parent(tree[offspring_key], key)
    return spouse_key == parent_key

def parent_is_sibling(tree, key, direction):
    '''
    key is an admixture node. This function checks if the parent in the direction of 'direction' is also the child of the parent in the direction of 
    'other_branch(direction)'. 
    '''
    parent_key=tree[key][direction]
    other_parent_key=tree[key][other_branch(direction)] #there is only one because key is an admixture node
    return (parent_key=='r' or other_parent_key in get_real_parents(tree[parent_key]))

    
    
        

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

def to_networkx_format(tree):
    edges=[]
    admixture_nodes=[]
    leaves=[]
    root=['r']
    coalescence_nodes=[]
    for key,node in tree.items():
        if node_is_leaf_node(node):
            leaves.append(key)
        else:
            if node_is_coalescence(node):
                coalescence_nodes.append(key)
            if node_is_admixture(node):
                admixture_nodes.append(key)
        ps=get_real_parents(node)
        for p in ps:
            edges.append((p,key))
    return leaves, admixture_nodes, coalescence_nodes, root, edges
    

def make_consistency_checks(tree, leaf_nodes=None):
    key_to_parents_by_def={key:[] for key in tree.keys()}
    key_to_children_by_def={key:[] for key in tree.keys()}
    key_to_children_by_family={key:[] for key in tree.keys()}
    key_to_parents_by_family={key:[] for key in tree.keys()}
    rooted_nodes=[]
    doppel_bands=[]
    pseudo_nodes=[]
    child_is_parent=[]
    recorded_leaf_nodes=[]
    illegal_admixture_props=[]
    illegal_branch_lengths=[]
    for key,node in tree.items():
        parents=[r for r in get_parents(node) if r is not None]
        children=[r for r in get_children(node) if r is not None]
        key_to_parents_by_def[key]+=[r for r in parents if r!='r']
        key_to_children_by_def[key]+=children
        if len(children)==2 and children[0]==children[1]:
            doppel_bands.append((key, ('children', children[0],children[1])))
        if len(parents)==2 and parents[0]==parents[1]:
            doppel_bands.append((key, ('parents', parents[0],parents[1])))
        if len(children)==1 and len(parents)==1:
            pseudo_nodes.append((key, ('child', children[0]), ('parent', parents[0])))
        if len(children)==0 and len(parents)==1:
            recorded_leaf_nodes.append(key)
        for r in children:
            key_to_parents_by_family[r]+=[key]
        for r in parents:
            if r!='r':
                key_to_children_by_family[r]+=[key]
            else:
                rooted_nodes.append(key)   
        if list(set(parents).intersection(children)):
            child_is_parent.append((key, ('parents', str(parents)), ('children', str(children))))
        if node_is_non_admixture(node):
            if node[3]<0:
                illegal_branch_lengths.append((key, (node[0], node[3])))
        else:
            if node[3]<0:
                illegal_branch_lengths.append((key, (node[0], node[3])))
            if node[4]<0:
                illegal_branch_lengths.append((key, (node[1], node[4])))
            if node[2]<0 or node[2]>1:
                illegal_admixture_props.append((key, node[2]))
                  
    
    def _transform_dic(dic):
        for key in dic.keys():
            dic[key]=set(dic[key])
        return dic
    key_to_children_by_def=_transform_dic(key_to_children_by_def)
    key_to_parents_by_def=_transform_dic(key_to_parents_by_def)
    key_to_children_by_family=_transform_dic(key_to_children_by_family)
    key_to_parents_by_family=_transform_dic(key_to_parents_by_family)
    def _print_inconsistencies(dic_def,dic_fam, pref=''):
        res=''
        for key in dic_def.keys():
            if dic_def[key]!=dic_fam[key]:
                res+=key+' '+pref +': '+str(dic_def[key])+'><'+str(dic_fam[key])+'  '
        return res
    
    family1=_print_inconsistencies(key_to_children_by_def, key_to_children_by_family, 'ch')
    family2=_print_inconsistencies(key_to_parents_by_def, key_to_parents_by_family, 'pa')
    
    bools=[]
    names=[]
    messages=[]
    
    
    consensus_bool=(key_to_children_by_def == key_to_children_by_family and key_to_parents_by_def == key_to_parents_by_family)
    consensus_message=family1+family2
    bools.append(consensus_bool)
    names.append('consensus')
    messages.append(consensus_message)
    
    roots_bool=(len(rooted_nodes)==2)
    roots_message=str(rooted_nodes)
    bools.append(roots_bool)
    names.append('roots')
    messages.append(roots_message)
    
    doppel_bands_bool=(len(doppel_bands)==0)
    doppel_bands_message=str(doppel_bands)
    bools.append(doppel_bands_bool)
    names.append('doppel_bands')
    messages.append(doppel_bands_message)
    
    pseudo_nodes_bool=(len(pseudo_nodes)==0)
    pseudo_nodes_message=str(pseudo_nodes)
    bools.append(pseudo_nodes_bool)
    names.append('pseudo_nodes')
    messages.append(pseudo_nodes_message)
    
    child_is_parent_bool=(len(child_is_parent)==0)
    child_is_parent_message=str(child_is_parent)
    bools.append(child_is_parent_bool)
    names.append('child_is_parent')
    messages.append(child_is_parent_message)
 
    if leaf_nodes is not None:
        leaf_nodes_bool=(set(leaf_nodes)==set(recorded_leaf_nodes))
    else:
        leaf_nodes_bool=True
    if leaf_nodes_bool:
        leaf_nodes_message=str(recorded_leaf_nodes)
    else:
        sl=set(leaf_nodes)
        srl=set(recorded_leaf_nodes)
        leaf_nodes_message=str(set(sl)-set(srl))+'><'+str(set(srl)-set(sl))
    bools.append(leaf_nodes_bool)
    names.append('leaf_nodes')
    messages.append(leaf_nodes_message) 
    
    illegal_branch_lengths_bool=(len(illegal_branch_lengths)==0)
    illegal_branch_lengths_message=str(illegal_branch_lengths)
    bools.append(illegal_branch_lengths_bool)
    names.append('illegal_branch_lengths')
    messages.append(illegal_branch_lengths_message)
    
    illegal_admixture_props_bool=(len(illegal_admixture_props)==0)
    illegal_admixture_props_message=str(illegal_admixture_props)
    bools.append(illegal_admixture_props_bool)
    names.append('illegal_admixture_props')
    messages.append(illegal_admixture_props_message)

    
    res_bool=all(bools)
    res_dic={name:(bool, message) for name,bool,message in zip(names,bools, messages)}
    
    return res_bool, res_dic



                
        
    


    
if __name__=='__main__':
        tree={'s1':['s1s2',None, None, 0.1,None],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3':['r',None, None, 0.4, None]}
        
        print create_burled_leaved_tree(5,1.0)
        print pretty_string(create_burled_leaved_tree(5,1.0))
        print pretty_string(create_balanced_tree(15, 1.0))
        
        from copy import deepcopy
        from Rcatalogue_of_trees import *
        
        tree2=insert_children_in_tree(tree_on_the_border2)
        print tree2
        print tree_on_the_border2_with_children
        
        #tree3=remove_parent_attachment(deepcopy(tree2), "s1")[0]
        #print tree3
        #print graft(tree3, 'f', 'r', 0.3, 'new_code', 0)
        #print remove_parent_attachment(tree2, 'f')
        #print get_descendants_and_rest(tree2, 'b')
        #print is_root(*get_parents(tree2['a']))
        
#         trouble={'a': ['132', 'c', 0.5, 0.06013128348912011, 0.1, 's2', None], 'c': ['212', 'r', 0.5, 0.03639623632910125, 0.15000000000000002, 'a', None], 'b': ['r', None, None, 0.07, None, '132', 's4'], 'e': ['132', None, None, 0.05, None, '212', 's3'], 's3': ['e', None, None, 0.3, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['212', None, None, 0.1, None, None, None], 's4': ['b', None, None, 0.3, None, None, None], '132': ['b', None, None, 0.1398687165108799, None, 'a', 'e'], '212': ['e', None, None, 0.06360376367089875, None, 'c', 's1']}
#         trouble= remove_parent_attachment(trouble, 'b')[0]
#         print trouble
#         print graft(trouble, 'b', 'r', 0.314, 'dont see me', 'hallo')
#         
#         trouble3={'a': ['n17', 'c', 0.5, 0.0006670327290825764, 0.1, 's2', None], 'c': ['n15', 'r', 0.5, 0.02087163982263861, 0.4814480657456043, 'a', None], 'n16': ['n17', None, None, 0.005272434567465561, None, 's4', 's3'], 'n17': [None, None, None, 0.013899593800954894, None, 'a', 'n16'], 'n15': ['r', None, None, 0.05969046586907494, None, 'c', 's1'], 's3': ['n16', None, None, 0.07815645814883887, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n15', None, None, 0.5947563021746359, None, None, None], 's4': ['n16', None, None, 0.00017898147838901196, None, None, None]}
#         print graft(trouble3, 'n17', 'a', 1, 'n18', 1)
#         
#         tree_trouble={'a': ['n37', 'c', 0.5, 1.5717637721311875, 0.1, 's1', None], 'n66': ['r', None, None, 0.008798782728668674, None, 's3', 's4'], 'c': ['n54', 'n37', 0.5, 0.771318479326775, 0.07345113788460944, 'a', None], 's3': ['n66', None, None, 0.010969920361510089, None, None, None], 's2': ['n54', None, None, 0.404441491678861, None, None, None], 's1': ['a', None, None, 0.06451508173696463, None, None, None], 's4': ['n66', None, None, 1.7305330689019498, None, None, None], 'n67': ['r', None, None, 0.24519067463109384, None, 'n54', 'n37'], 'n54': ['n67', None, None, 0.25870104556004564, None, 'c', 's2'], 'n37': ['n67', None, None, 0.9342460567572629, None, 'c', 'a']}
#         print 'tree_trouble', tree_trouble
#         removed=remove_parent_attachment(tree_trouble, 'n54')[0]
#         print 'pruned tree', removed
#         adm=graft(removed, 'n54', 'c', 0.3, 'n68', 0) #FIXME: the function is being called like this via make_regraft
#         print adm
#         
#         print tree_on_the_border2_with_children
#         print remove_admix(tree_on_the_border2_with_children, 'a', 0)
#         
#         trouble_tree={'157': ['e', '95', 0.48, 0.053593995338132035, 0, '209', None], '156': ['r', None, None, 0.03852924467461515, None, '177', '184'], '196': ['209', '184', 0.48, 0.0005592709077173309, 0, 'c', None], '177': ['156', '87', 0.48, 0.04409469700439552, 0, '226', None], '138': ['54', '184', 0.48, 0.0923842952243634, 0, 's3', None], '87': ['95', None, None, 0.021717663805443786, None, 'b', '177'], '251': ['c', 'r', 0.48, 0.048451187058856385, 0, 'a', None], 's3': ['138', None, None, 0.003662331553773855, None, None, None], 's2': ['149', None, None, 0.031232370021655256, None, None, None], 's1': ['d', None, None, 0.1, None, None, None], 's4': ['b', None, None, 0.3, None, None, None], '184': ['156', None, None, 0.014531521069126824, None, 'f', '196'], '209': ['157', '164', 0.48, 0.04266174601860901, 0, '196', None], '149': ['a', None, None, 0.018767629978344746, None, 's2', '54'], '164': ['f', None, None, 0.00649759823036904, None, '95', '209'], '226': ['177', None, 0.48, 0.004168053034161371, 0, 'd', None], '95': ['164', None, None, 0.02046617829442733, None, '87', '157'], 'a': ['b', '251', 0.5, 0.2, 0.05154881294114362, '149', None], 'c': ['196', 'd', 0.5, 0.003184987735541627, 0.1, '251', None], 'b': ['87', None, None, 0.001318559669759849, None, 's4', 'a'], 'e': ['f', None, None, 0.05, None, '157', '184'], 'd': ['226', None, None, 0.0017372499614431128, None, 's1', 'c'], 'f': ['184', None, None, 0.005468478930873175, None, '164', 'e'], '54': ['184', '149', 0.48, 0.17449964855898192, 0, '138', None]}

        
        
        print make_consistency_checks(tree_good)
        print make_consistency_checks(tree_without_consensus)
        print make_consistency_checks(tree_with_self_connection)
        print make_consistency_checks(tree_with_pseudo_node)
        print make_consistency_checks(tree_with_doppel_band)
        print make_consistency_checks(tree_with_negative_bl)
        print make_consistency_checks(tree_with_illegal_alpha)
        
        #print update_all_branches(tree_good, )
        
        from piece_of_tree import piece
        p=piece(0, None, 0, None, 'r', 0, None)
        pretty_print(  move_node(tree_good, 'c', 1, 'd', 0.01, p, new_node_name='x')  )

        
        
        