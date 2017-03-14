
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
    '''
    constructs tree of the form 
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
    root_keys=_find_rooted_nodes(tree)

    if len(root_keys)==1:#this is the special case where an admixture leads to the new root
        return graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_keys[0], remove_branch=remove_branch)
    
    #updating the grafted_branch. Easy.
    tree[remove_key][remove_branch]='r'
    
    #dealing with the other child of the new node, but since the new node is the root, the old root is the new node. If that makes sense.
    #print 'root_keys', root_keys
    tree[new_name_for_old_root]=['r', None, None, insertion_spot, None,root_keys[0][0], root_keys[1][0]]
    
    #dealing with the children of the new node.
    for rkey, b_l in root_keys:
        tree[rkey]=_update_parent(tree[rkey], 'r', new_name_for_old_root)
    
    return tree

def graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_key, remove_branch=0):
    print 'undoing a closed branch', insertion_spot, remove_key, root_key
    tree[remove_key][remove_branch]='r'
    tree[root_key[0]],_=get_branch_length_and_reset(tree[root_key[0]], 'r', insertion_spot)
    return tree

def get_number_of_admixes(tree):
    return sum((1 for node in tree.values() if node_is_admixture(node)))
            
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
    if parent_key!='r':
        tree[parent_key]=_update_child(tree[parent_key], old_child=rkey, new_child=orphanota_key)
    if source_key=='r':
        tree,_,t3=remove_root_attachment(tree, rkey) #now sorphanota_key is new root
        del tree[rkey]
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
            edges.append((key,p))
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
    
    
    
    
    res_bool=all(bools)
    res_dic={name:(bool, message) for name,bool,message in zip(names,bools, messages)}
    
    return res_bool, res_dic

                
        
    


    
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
        
        print tree_on_the_border2_with_children
        print remove_admix(tree_on_the_border2_with_children, 'a', 0)
        
        trouble_tree={'157': ['e', '95', 0.48, 0.053593995338132035, 0, '209', None], '156': ['r', None, None, 0.03852924467461515, None, '177', '184'], '196': ['209', '184', 0.48, 0.0005592709077173309, 0, 'c', None], '177': ['156', '87', 0.48, 0.04409469700439552, 0, '226', None], '138': ['54', '184', 0.48, 0.0923842952243634, 0, 's3', None], '87': ['95', None, None, 0.021717663805443786, None, 'b', '177'], '251': ['c', 'r', 0.48, 0.048451187058856385, 0, 'a', None], 's3': ['138', None, None, 0.003662331553773855, None, None, None], 's2': ['149', None, None, 0.031232370021655256, None, None, None], 's1': ['d', None, None, 0.1, None, None, None], 's4': ['b', None, None, 0.3, None, None, None], '184': ['156', None, None, 0.014531521069126824, None, 'f', '196'], '209': ['157', '164', 0.48, 0.04266174601860901, 0, '196', None], '149': ['a', None, None, 0.018767629978344746, None, 's2', '54'], '164': ['f', None, None, 0.00649759823036904, None, '95', '209'], '226': ['177', None, 0.48, 0.004168053034161371, 0, 'd', None], '95': ['164', None, None, 0.02046617829442733, None, '87', '157'], 'a': ['b', '251', 0.5, 0.2, 0.05154881294114362, '149', None], 'c': ['196', 'd', 0.5, 0.003184987735541627, 0.1, '251', None], 'b': ['87', None, None, 0.001318559669759849, None, 's4', 'a'], 'e': ['f', None, None, 0.05, None, '157', '184'], 'd': ['226', None, None, 0.0017372499614431128, None, 's1', 'c'], 'f': ['184', None, None, 0.005468478930873175, None, '164', 'e'], '54': ['184', '149', 0.48, 0.17449964855898192, 0, '138', None]}

        tree_good={'s1':['d',None, None, 0.1,None,None,None],
              's2':['a',None, None,0.05,None,None,None],
              's3':['e',None,None, 0.3,None,None,None],
              's4':['b',None,None, 0.3,None,None,None],
              'a':['b','c', 0.5,0.2,0.1,'s2',None],
              'c':['e','d',0.5,0.1,0.1,'a',None],
              'b':['f',None,None,0.05,None,'s4','a'],
              'f':['r',None,None,0.02,None,'b','e'],
              'e':['f',None,None,0.05,None,'c','s3'],
              'd':['r',None,None,0.05,None,'s1','c']}

        tree_without_consensus={'s1':['d',None, None, 0.1,None,None,None],
              's2':['a',None, None,0.05,None,None,None],
              's3':['e',None,None, 0.3,None,None,None],
              's4':['a',None,None, 0.3,None,None,None],
              'a':['b','c', 0.5,0.2,0.1,'s2',None],
              'c':['e','d',0.5,0.1,0.1,'a',None],
              'b':['f',None,None,0.05,None,'s4','a'],
              'f':['r',None,None,0.02,None,'b','e'],
              'e':['f',None,None,0.05,None,'c','s3'],
              'd':['r',None,None,0.05,None,'s1','c']}
        
        tree_with_self_connection={'s1':['r',None, None, 0.1,None,None,None],
              's2':['a',None, None,0.05,None,None,None],
              's3':['f',None,None, 0.3,None,None,None],
              's4':['b',None,None, 0.3,None,None,None],
              'a':['b','c', 0.5,0.2,0.1,'s2',None],
              'c':['c',None, 0.5,0.1,0.1,'a','c'],
              'b':['f',None,None,0.05,None,'s4','a'],
              'f':['r',None,None,0.02,None,'b','s3']}
        
        tree_with_pseudo_node={'s1':['d',None, None, 0.1,None,None,None],
              's2':['a',None, None,0.05,None,None,None],
              's3':['e',None,None, 0.3,None,None,None],
              's4':['b',None,None, 0.3,None,None,None],
              'a':['b','c', 0.5,0.2,0.1,'s2',None],
              'c':['e',None,None,0.1,None,'a',None],
              'b':['f',None,None,0.05,None,'s4','a'],
              'f':['r',None,None,0.02,None,'b','e'],
              'e':['f',None,None,0.05,None,'c','s3'],
              'd':['r',None,None,0.05,None,'s1',None]}
        
        tree_with_doppel_band={'s1':['d',None, None, 0.1,None,None,None],
              's2':['a',None, None,0.05,None,None,None],
              's3':['d',None,None, 0.3,None,None,None],
              's4':['b',None,None, 0.3,None,None,None],
              'a':['b','c', 0.5,0.2,0.1,'s2',None],
              'c':['e','e',0.5,0.1,0.1,'a',None],
              'b':['f',None,None,0.05,None,'s4','a'],
              'f':['r',None,None,0.02,None,'b','e'],
              'e':['f',None,None,0.05,None,'c','c'],
              'd':['r',None,None,0.05,None,'s1','s3']}
        
        print make_consistency_checks(tree_good)
        print make_consistency_checks(tree_without_consensus)
        print make_consistency_checks(tree_with_self_connection)
        print make_consistency_checks(tree_with_pseudo_node)
        print make_consistency_checks(tree_with_doppel_band)

        
        
        