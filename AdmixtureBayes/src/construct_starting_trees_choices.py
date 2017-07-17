from generate_prior_trees import generate_phylogeny, simulate_number_of_admixture_events
from tree_statistics import identifier_to_tree_clean
from Rtree_operations import create_trivial_tree

from copy import deepcopy
from warnings import warn

def get_starting_trees(inputs, no_chains, random=False, nodes=None, skewed_admixture_prior=False):
    
    if isinstance(inputs, basestring):
        tree=input_to_tree(inputs, random, nodes, skewed_admixture_prior)
        return [deepcopy(tree) for _ in xrange(no_chains)]
    else:
        list_of_trees=[input_to_tree(input, random, nodes, skewed_admixture_prior) for input in inputs]
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
                return identifier_to_tree_clean(f.readline().rstrip(), nodes=nodes)
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
                stree=create_trivial_tree(n, tree)
    else:#is it a tree already?
        return input