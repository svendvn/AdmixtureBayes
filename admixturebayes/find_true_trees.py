from tree_statistics import identifier_to_tree_clean, unique_identifier_and_branch_lengths
from post_analysis import read_tree_file
from Rtree_operations import change_admixture, get_categories, get_leaf_keys
from copy import deepcopy

class tree_unifier(object):
    
    def __init__(self):
        '''
        key:(val1,val2, val3) where 
            key=lookup topology string, 
            val1 is unique topology string, 
            val2 is the permutation of branches and 
            val3 is the (signed) permutation of the admixture proportions.
        '''
        self.seen_trees={}
        
    def __call__(self, stree):
        topology,branches,admixtures=stree.split(';')
        if topology in self.seen_trees:
            target_topology, branch_permutation, admixture_permutation= self.seen_trees[topology]
        else:
            update_dic= analyze_tree(topology, branches, admixtures)
            if len(update_dic)>100: #computational reasons
                self.seen_trees[topology]=update_dic[topology]
            else:
                self.seen_trees.update(update_dic)
            #print self.seen_trees
            target_topology, branch_permutation, admixture_permutation=self.seen_trees[topology]
        #print len(self.seen_trees)
        new_branch_string=make_branch_string(branches, branch_permutation)
        new_admixtures_string=make_admixture_string(admixtures, admixture_permutation)
        return ';'.join([target_topology, new_branch_string, new_admixtures_string])
    
def make_branch_string(branches, branch_permutations):
    branch_pieces=branches.split('-')
    return '-'.join([branch_pieces[branch_permutations[i]] for i in range(len(branch_pieces))])

def make_admixture_string(admixes, admixture_permutations):
    if not admixes:
        return ''
    res=[]
    admixture_pieces=admixes.split('-')
    for i in range(len(admixture_pieces)):
        target_admixture=admixture_permutations[i]
        if target_admixture<-0.5:
            res.append(1.0-float(admixture_pieces[abs(target_admixture)-1]))
        else:
            res.append(float(admixture_pieces[abs(target_admixture)-1]))
    return '-'.join(map(str,res))
        

def make_possible_files(true_tree_file, res_file):
    tree, nodes=read_tree_file(true_tree_file)
    possible_strees=get_possible_strees(tree, nodes)
    with open(res_file, 'w') as f:
        f.write(' '.join(nodes)+'\n')
        for stree in possible_strees:
            f.write(stree+'\n')
            
def get_unique_plottable_tree(tree, nodes=None):
    if nodes is None:
        nodes=sorted(get_leaf_keys(tree))
    possible_strees=sorted(get_possible_strees(tree, nodes))
    return possible_strees[0]
    
def get_possible_strees(tree, nodes):
    
    leaves,_,admixture_keys=get_categories(tree)
    k=len(admixture_keys)
    format_code='{0:0'+str(k)+'b}'
    
    
    n_trees=[]
    for i in range(2**k):
        pruned_tree = deepcopy(tree)
        bina= format_code.format(i)
        prop=1.0
        for adm_key,str_bin in zip(admixture_keys, list(bina)):
            int_bin=int(str_bin)
            if int_bin==1:
                pruned_tree[adm_key]=change_admixture(pruned_tree[adm_key])
        n_tree= unique_identifier_and_branch_lengths(pruned_tree, leaf_order=nodes)
        n_trees.append(n_tree)
    return n_trees


def permutation_stree_to_permutations(permutation_stree):
    topology, branches, admixtures=permutation_stree.split(';')
    branches_permutation=invert(branches.split('-'))
    admixtures_signs=[-1 if i<0 else 1 for i in map(int, admixtures.split('-'))]
    admixtures_values=map(abs, map(int, admixtures.split('-')))
    admixtures_permutation=invert(admixture_values)

def analyze_tree(topology, branches, admixtures):
    
    id_branches='-'.join(map(str, range(len(branches.split('-')))))
    id_admixtures='-'.join(map(str, range(1,len(admixtures.split('-'))+1)))
    #print branches, id_branches
    id_stree=';'.join([topology,id_branches, id_admixtures])
    no_leaves=len((id_stree.split('-')[0]).split('.'))
    id_tree=identifier_to_tree_clean(id_stree.strip())
    
    strees= sorted(get_possible_permutation_strees(id_tree))
    top_topology=strees[0].split(';')[0]
    res={}
    for stree in strees:
        lookup_topology, branches_sperm, admixtures_sperm= stree.split(';')
        rf=map(round, map(float, branches_sperm.split('-')))
        branches_permutation= map(int, rf)
        admixtures_permutation=get_admixtures_permutation(admixtures_sperm)
        res[lookup_topology]=(top_topology, branches_permutation, admixtures_permutation)
    #print res
    return res
    
def get_admixtures_permutation(admixtures):
    res=[]
    for a in admixtures.split('-'):
        if a:
            number=int(round(float(a)))
            if number>100000:
                res.append(-(number%100000))
            else:
                res.append(number)
    return res
    

def get_possible_permutation_strees(tree):
    leaves,_,admixture_keys=get_categories(tree)
    k=len(admixture_keys)
    format_code='{0:0'+str(k)+'b}'
    
    n_trees=[]
    for i in range(2**k):
        pruned_tree = deepcopy(tree)
        bina= format_code.format(i)
        prop=1.0
        for adm_key,str_bin in zip(admixture_keys, list(bina)):
            int_bin=int(str_bin)
            if int_bin==1:
                pruned_tree[adm_key]=change_admixture(pruned_tree[adm_key])
        n_tree= unique_identifier_and_branch_lengths(pruned_tree)
        n_trees.append(n_tree)
    return n_trees
    

def change_admixture_permutation(node):
    new_node=[node[1],node[0],100000+node[2],node[4],node[3]]+node[5:]
    return new_node

if __name__=='__main__':
    tu=tree_unifier()
    
    print tu('a.w.w-c.0.c.2-c.0;1-1-1-1-1-1-1;0.3')
    print tu('a.w.w-c.c.0.1-c.0;1-1-1-1-1-1-1;0.3')
    
    import sys
    sys.exit()
    strees=[]
    with open('example_strees.txt', 'r') as f:
        for l in f.readlines():
            strees.append(l.strip())
    for stree in strees:
        print stree, tu(stree)
    import sys
    sys.exit()
    from argparse import ArgumentParser

    parser = ArgumentParser(usage='pipeline for Admixturebayes', version='1.0.0')

    #overall options
    parser.add_argument('--i', type=str, default='_true_tree.txt', help='')
    parser.add_argument('--o', type=str, default='tree.txt')
    
    options=parser.parse_args()
    make_possible_files(options.i, options.o)