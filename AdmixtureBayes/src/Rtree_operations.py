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

tree_admix_to_child={
    's1':['r',None,None, 0.1,None],
    's2':['s2a',None,None,0.1,None],
    's3':['s3s2',None,None,0.1,None],
    's2a':['s3s2','s3s2a', 0.5,0.1,0.13],
    's3s2':['s3s2a',None,None,0.1,None],
    's3s2a':['r',None,None,0.01]
    }

def regraft(tree, remove_key, add_to_branch):
    tree=remove_parent_attachment(tree, remove_key)
    
def remove_parent_attachment(tree, orphanonte_key):
    '''
    This takes the tree and removes the parent of orphanonte_key.
    '''
    pkey=get_parents(tree[orphanonte_key])[0]
    children_of_parent=get_children(tree[pkey])
    for key in children_of_parent:
        if key is not None:
            tree[key]=
    

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

def get_children(node):
    return node[5:7]

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
    children={key:[None, None] for key in tree}
    for key in tree:
        parent1, parent2 = get_parents(tree[key])
        print parent1, parent2, key
        print children
        if parent1!='r':
            if children[parent1][0] is None:
                children[parent1][0] = key
            else:
                children[parent1][1]=key
        if parent2 is not None and parent2!='r':
            if children[parent2][0] is None:
                children[parent2][0] = key
            else:
                children[parent2][1]=key
    for key in tree:
        tree[key]=tree[key][:5]+children[key]
    return tree
        
def get_parents(node):
    return node[:2]

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
        
        print get_children_and_rest(tree2, 'b')
        print is_root(*get_parents(tree2['a']))
        