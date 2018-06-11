from argparse import ArgumentParser
from tree_statistics import identifier_file_to_tree_clean, unique_identifier_and_branch_lengths
from Rtree_operations import get_leaf_keys, add_outgroup, pretty_string
from subgraphing import get_subtree
from tree_to_data import unzip, gzip
import pandas as pd
from copy import deepcopy



usage = """%(prog)s [options] <forwarder dirs>

This program estimates the coalescence and migration rates over time together with a constant
recombination rate."""

parser = ArgumentParser(usage=usage, version="%(prog)s 1.0")

parser.add_argument("--input_file",
                    type=str,
                    default="",
                    help="")

parser.add_argument('--input_add', type=str, default='', help='If tree options is also chosen as input_type this should be a file with distance to the outgroup. Otherwise assumed to be 0')
parser.add_argument('--outgroup_name', type=str, default='out', help='If input_add is supplied you would also need to specify the name of the outgroup. That is done with this argument.')

parser.add_argument('--input_type', type=str, default='snps', choices=['snps','tree'], help='The program can make both snp files and tree files smaller. This tells which type')

parser.add_argument('--populations', type=str, nargs='+', help='The populations to keep after having removed the other. If a population is not found it will raise an error.')
parser.add_argument('--output_file', type=str, default='', help='The file where the populations should be saved.')

options= parser.parse_args()

if options.input_type=='tree':
    tree=identifier_file_to_tree_clean(options.input_file)
    if options.input_add:
        with open(options.input_add, 'r') as f:
            add=float(f.readline())
        tree=          add_outgroup(tree,  
                                    inner_node_name='new_node', 
                                    to_new_root_length=float(add),
                                    to_outgroup_length=0, 
                                    outgroup_name=options.outgroup_name)
    nodes=get_leaf_keys(tree)
    assert all((a in nodes for a in options.populations)), 'Requested population was not found in the tree'
    subtree=get_subtree(tree, options.populations)
    if not options.output_file:
        options.output_file=options.input_file+'_'.join(options.populations)
    with open(options.output_file, 'w') as f:
        f.write(' '.join(sorted(options.populations))+'\n')
        f.write(unique_identifier_and_branch_lengths(subtree))
if options.input_type=='snps':
    if options.input_file.endswith('.gz'):
        options.input_file=unzip(options.input_file, overwrite=True)
    df=pd.read_csv(options.input_file, usecols=options.populations, sep=' ')
    if not options.output_file:
        options.output_file=options.input_file+'_'.join(options.populations)
    df.to_csv(options.output_file, sep=' ', index=False)
    gzip(options.output_file, overwrite=True)



    
    
    