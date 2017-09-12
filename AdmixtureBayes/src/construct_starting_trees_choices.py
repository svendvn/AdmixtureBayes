from generate_prior_trees import generate_phylogeny, simulate_number_of_admixture_events
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from Rtree_operations import create_trivial_tree, scale_tree

from copy import deepcopy
from warnings import warn

def get_starting_trees(inputs, no_chains, random=False, nodes=None, skewed_admixture_prior=False, add=0, scaling=0):
    if add:
        with open(add, 'r') as f:
            add=float(f.readline())
    else:
        add=0
    if isinstance(inputs, basestring):
        tree=input_to_tree(inputs, random, nodes, skewed_admixture_prior)
        if scaling>1e-9:
            scale_tree(tree, scaling)
            addt=add*scaling
        return [(deepcopy(tree),addt) for _ in xrange(no_chains)]
    else:
        list_of_trees=[]
        for input in inputs:
            tree=input_to_tree(input, random, nodes, skewed_admixture_prior)
            if scaling>1e-9:
                scale_tree(tree, scaling)
                addt=add*scaling
            list_of_trees.append((tree,addt))
        if len(list_of_trees)>no_chains:
            warn('Some starting trees will not be used')
        starting_trees=[]
        count=0
        for _ in xrange(no_chains):
            starting_trees.append(list_of_trees[count])
            count+=1
            if count==len(list_of_trees):
                count=0
    return starting_trees
            
    
    
def input_to_tree(input, random=False, nodes=None, skewed_admixture_prior=False):
    if isinstance(input, basestring):
        if ';' in input:
            return identifier_to_tree_clean(input, nodes=nodes)
        elif '.' in input:
            with open(input, 'r') as f:
                f.readline() #removing empty file
                return identifier_to_tree_clean(f.readline().rstrip(), leaves=generate_predefined_list_string(deepcopy(nodes)))
        elif ',' in input:
            first_part, last_part= input.split(',')
            n=int(first_part[1:])
            k=int(last_part[:-1])
            return generate_phylogeny(n,k, nodes=nodes, skewed_admixture_prior=skewed_admixture_prior)
        else:
            n=int(input)
            if random:
                return generate_phylogeny(n,nodes=nodes, skewed_admixture_prior=skewed_admixture_prior)
            else:
                return create_trivial_tree(n, 1.0)
    else:#is it a tree already?
        return input