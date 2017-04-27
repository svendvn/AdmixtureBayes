from Rtree_operations import find_rooted_nodes, get_real_children, get_real_parents, get_branch_length, mother_or_father


class piece(object):
    
    def __init__(self, start_lattitude, end_lattitude, start_distance, end_distance, child_key, child_branch, parent_key):
        self.start_lattitude=start_lattitude
        self.end_lattitude=end_lattitude
        self.start_distance=start_distance
        self.end_distance=end_distance
        self.child_key=child_key
        self.child_branch=child_branch
        self.parent_key=parent_key
        if start_lattitude>end_lattitude:
            self.direction='to_leaves'
        else:
            self.direction='to_root'
        
    def __str__(self):
        return ', '.join(map(str,[self.start_lattitude, self.end_lattitude, self.start_distance, self.end_distance]))

class lineage(object):
    
    def __init__(self, key,  distance=0, lattitude=0):
        self.key,self.distance, self.lattitude= key, distance, lattitude
        
    def follow(self, tree, visited_keys=[]):
        node=tree[self.key]
        new_lineages=[]
        pieces=[]
        new_keys_visited=[]
        for key in get_real_children(tree[self.key]):
            if key not in visited_keys:
                branch=mother_or_father(tree, child_key=key, parent_key=self.key)
                l=get_branch_length(tree, key, branch)
                pieces.append(piece(self.lattitude, self.lattitude-l, self.distance, self.distance+l, key, branch, self.key))
                new_lineages.append(lineage(key,self.distance+l, self.lattitude-l))
        for key in get_real_parents(tree[self.key]):
            if key not in visited_keys:
                branch=mother_or_father(tree, child_key=self.key, parent_key=key)
                l=get_branch_length(tree, self.key, branch)
                pieces.append(piece(self.lattitude, self.lattitude+l, self.distance, self.distance+l, self.key, branch, key))
                new_lineages.append(lineage(key,self.distance+l, self.lattitude+l))
        if self.key=='r':#add the very long piece
            pieces.append(piece(self.lattitude, None, self.distance, None,'r',0,None))
        return new_lineages, pieces


def insert_root_in_tree(tree):
    (child_key1,_,_),(child_key2,_,_)=find_rooted_nodes(tree)
    tree['r']=[None,None,None,None,None,child_key1, child_key2]

def remove_root_from_tree(tree):
    del tree['r']
    

def distanced_branch_lengths(tree, start_key):
    insert_root_in_tree(tree)
    pieces=[]
    visited_keys=[]
    lineages=[lineage(start_key, 0, 0)]
    while lineages:
        lineages.sort(key=lambda x: x.distance)
        lin=lineages.pop(0)
        if lin.key not in visited_keys:
            visited_keys.append(lin.key)
            print lin.key
            new_lineages, new_pieces= lin.follow(tree, visited_keys)
            pieces.extend(new_pieces)
            lineages.extend(new_lineages)
    remove_root_from_tree(tree)
    return pieces
    
if __name__=='__main__':
    from Rcatalogue_of_trees import tree_good
    ad=distanced_branch_lengths(tree_good, 'a')
    for e in ad:
        print e
    