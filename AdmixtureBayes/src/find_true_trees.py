from tree_statistics import identifier_to_tree_clean, unique_identifier_and_branch_lengths
from post_analysis import read_tree_file
from Rtree_operations import change_admixture, get_categories
from copy import deepcopy


def make_possible_files(true_tree_file, res_file):
    tree, nodes=read_tree_file(true_tree_file)
    possible_strees=get_possible_strees(tree, nodes)
    with open(res_file, 'w') as f:
        f.write(' '.join(nodes)+'\n')
        for stree in possible_strees:
            f.write(stree+'\n')
    
        
    
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

if __name__=='__main__':
    from tree_generation_laboratory import simulate_tree, save_tree_and_nodes
    from Rtree_operations import get_trivial_nodes
    stree=simulate_tree(6,2)
    save_tree_and_nodes(stree, 'tmp.txt', get_trivial_nodes(6))
    make_possible_files('tmp.txt', 'tmp_out.txt')