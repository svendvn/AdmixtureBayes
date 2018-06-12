from Rtree_to_covariance_matrix import get_populations
from Treemix_to_AdmixtureBayes import Node
from construct_nodes_choices import read_one_line
import pandas as pd
from argparse import ArgumentParser
from tree_statistics import identifier_to_tree_clean,generate_predefined_list_string,topological_identifier_to_tree_clean
from copy import deepcopy
#from tree_plotting import plot_node_structure_as_directed_graph, plot_as_directed_graph NOTICE THAT THIS IS CALLED ELSEWHERE!! IN THE SCRIPT
import sys

parser = ArgumentParser(usage='pipeline for consensus tree maker', version='1.0.0')

parser.add_argument('--test_run', default=False, action='store_true', help='will overwrite everything and run a test function')
parser.add_argument('--input_file', default='result_mc3.csv', type=str, help='The input file that should contain a column named tree or the option no_header should be turned on in which case every line is assumed to hold just one tree')
parser.add_argument('--nodes', default='', type=str, help='file where the first line is the leaf nodes')
parser.add_argument('--no_sort', default=False, action='store_true', help='often the tree is sorted according to the leaf names. no_sort willl assumed that they are not sorted according to this but sorted according to ')
parser.add_argument('--no_header', default=False, action='store_true',help='will assume that there is no header in the file')
parser.add_argument('--burn_in_rows', default=0, type=int, help='the number of rows that will be skipped in the input file as burn-in period')
parser.add_argument('--burn_in_fraction', default=0.0, type=float, help='the proportion of the rows that are discarded as burn in period')
parser.add_argument('--tree_column_name', default='tree', type=str, help='the name in the header of the column with all the trees.')
parser.add_argument('--max_number_of_trees', default=10000, type=int, help='an upper limit on the number of trees to reduce computational pressure')
parser.add_argument('--sep', default=',', type=str, help='the separator used in the input file')
parser.add_argument('--consensus_method', choices=['descendant_frequencies'], default='descendant_frequencies', help='Which method should be used to calculate the consensus tree?')
parser.add_argument('--min_w', default=0.0, type=float, help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')

parser.add_argument('--posterior_threshold', default=[0.25,0.5,0.75,0.9,0.95,0.99], type=float, nargs='+', help='The posterior threshold at which to include ')

parser.add_argument('--plot_tops_file', action='store_true', default=False, help='this will assume that the file is a tops file from downstream_analysis_parser and plot each line numbered.')

parser.add_argument('--get_effective_number_of_admixtures', action='store_true', default=False, help='this will cancel all the other analysis and only print the topological number of admixes(tadmixes) to a a file.')
parser.add_argument('--effective_number_of_admixtures_file', type=str, default='no_tadmixes.txt', help='this is the file in which to write the effective number of admixes in the file')
parser.add_argument('--suppress_plot', default=False, action='store_true')
options= parser.parse_args()



def combine_nodes(node_structure, new_node, seen_sets):
    candidate=new_node.name
    seen=[]
    for lists_of_fixed_size in seen_sets[::-1]:
        for attached_branch in lists_of_fixed_size:
            if( attached_branch.issubset(candidate) and 
               ((not attached_branch.issubset(seen)) or (not node_structure[attached_branch].has_parent()))):
                seen.extend(list(attached_branch))
                new_node.add_child(node_structure[attached_branch])
                node_structure[attached_branch].add_parent(new_node)
    return node_structure
                
                
def get_number_of_tadmixtures(node_structure):
    total=0
    for key in node_structure:
        total+=max(0,node_structure[key].get_number_of_parents()-1)
    return total
    
def node_combinations_to_node_structure(node_combinations):
    length_sorted={}
    for node_combination in node_combinations:
        leaves=frozenset(node_combination.split('.'))
        k=len(leaves)
        if k in length_sorted:
            length_sorted[k].append(leaves)
        else:
            length_sorted[k]=[leaves]
    length_sorted_list=[length_sorted.get(k,[]) for k in range(1,max(length_sorted.keys())+1)]
    #length_sorted_list is of the form [[[A],[B],[C]],[[A,B],[B,C]],...,[[A,B,C]]]
    node_structure={}
    for leaf_node in length_sorted_list[0]:
        node_structure[leaf_node]=Node(leaf_node)
    added_sets=[length_sorted_list[0]]
    for lists_of_fixed_size in length_sorted_list[1:]:
        for branch_set in lists_of_fixed_size:
            new_node=Node(branch_set)
            combine_nodes(node_structure, new_node, added_sets)
            node_structure[branch_set]=new_node
        added_sets.append(lists_of_fixed_size)
    return node_structure
            
        
if options.plot_tops_file:
    with open(options.input_file, 'r') as f:
        for lin in f.readlines():
            rank, probability, combination=lin.rstrip().split(',')
            all_nodes=[c.split('.') for c in combination.split('_')]
            flattened=[item for sublist in all_nodes for item in sublist]
            a=list(set(flattened))
            code=rank+'_'+str(int(100*round(float(probability),2)))+'_'+'_'.join(a)
            print 'code',code
            node_structure=node_combinations_to_node_structure(combination.split('_'))
            print node_structure
            plot_node_structure_as_directed_graph(node_structure, drawing_name=code+'.png')
    sys.exit()
    

if options.test_run:
    from generate_prior_trees import generate_phylogeny
    from tree_statistics import unique_identifier_and_branch_lengths
    from tree_plotting import plot_node_structure_as_directed_graph, plot_as_directed_graph 
    N=5
    tree1=generate_phylogeny(N,1)
    plot_as_directed_graph(tree1, drawing_name='tree1.png')
    tree2=generate_phylogeny(N,1)
    plot_as_directed_graph(tree2, drawing_name='tree2.png')
    stree1=unique_identifier_and_branch_lengths(tree1)
    stree2=unique_identifier_and_branch_lengths(tree2)
    with open('tmp_tree.txt','w') as f:
        f.write(' '.join(['s'+str(i) for i in range(1,N+1)])+'\n')
        f.write(stree1)
    with open('trees.txt','w') as f:
        f.write(stree1+'\n'+stree2+'\n'+stree1)
    
    options.input_file='trees.txt'
    options.nodes='tmp_tree.txt'
    options.no_header=True
    options.posterior_threshold=[0.25,0.5,0.9]
    
    
print 'Reading file...'
#loading trees
if options.no_header:
    strees=[]
    with open(options.input_file, 'r') as f:
        for lin in f.readlines():
            strees.append(lin.rstrip())
else:
    df=pd.read_csv(options.input_file, sep=options.sep, usecols=[options.tree_column_name])
    strees=df[options.tree_column_name].tolist()
n=len(strees)
print 'trees read: ',n

#thinning tree list
    
rows_to_remove_from_fraction=int(options.burn_in_fraction*n)
rows_to_remove=max(rows_to_remove_from_fraction, options.burn_in_rows)
strees=strees[rows_to_remove:]

print 'removed burn-in:', rows_to_remove
print 'In list are now', len(strees),'trees'

#thinning

distance_between=max(1,len(strees)//options.max_number_of_trees)
nstrees=[]
for a,stree in enumerate(strees):
    if a%distance_between==0 and len(nstrees)<options.max_number_of_trees:
        nstrees.append(stree)
print 'thinned'
print 'In list are now', len(nstrees),'trees'

N=len(nstrees)

seen_node_combinations={}

nodes=read_one_line(options.nodes)
if not options.no_sort:
    nodes=sorted(nodes)

    
tenth=len(nstrees)//10
for i,stree in enumerate(nstrees):
    if tenth>0 and i%tenth==0:
        print i//tenth*10, '%'
    if ';' in stree:
        tree=identifier_to_tree_clean(stree, leaves=generate_predefined_list_string(deepcopy(nodes)))
    else:
        tree=topological_identifier_to_tree_clean(stree, leaves=generate_predefined_list_string(deepcopy(nodes)))
    ad=get_populations(tree, min_w=options.min_w)
    for a in ad:
        seen_node_combinations[a]=seen_node_combinations.get(a,0)+1

for threshold in options.posterior_threshold:
    total_threshold=int(N*threshold)
    final_node_combinations=[k for k,v in seen_node_combinations.items() if v > total_threshold]
    print 'final_node_combinations', final_node_combinations
    final_node_structure=node_combinations_to_node_structure(final_node_combinations)
    if options.get_effective_number_of_admixtures:
        with open(options.effective_number_of_admixtures_file, 'w') as f:
            effictive_admixtures=get_number_of_tadmixtures(final_node_structure)
            f.write(str(effictive_admixtures))
    if not options.suppress_plot:
        from tree_plotting import plot_node_structure_as_directed_graph, plot_as_directed_graph 
        plot_node_structure_as_directed_graph(final_node_structure, drawing_name='tmp'+str(total_threshold)+'.png')
    
    


if __name__=='__main__':
    from generate_prior_trees import generate_phylogeny
    
    


    

