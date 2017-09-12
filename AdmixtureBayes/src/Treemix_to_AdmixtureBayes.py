import subprocess
from newick import parse_tree
from Rtree_operations import pretty_string, insert_children_in_tree
from meta_proposal import new_node_naming_policy

class node(object):
    
    def __init__(self):
        self.count=0
        
    def __call__(self):
        self.count+=1
        print 'n'+str(self.count)
        return 'n'+str(self.count)

def read_treemix_file(filename):
    np=new_node_naming_policy()
    reduced_filename='.'.join(filename.split(".")[:-1])
    take_copy_args=['cp', filename, filename+".tmp"]
    move_back_args=['mv', filename+'.tmp', filename]
    args=['gunzip', '-f', filename]
    subprocess.call(take_copy_args)
    subprocess.call(args)
    subprocess.call(move_back_args)
    with open(reduced_filename, 'r') as f:
        newick_tree=f.readline().rstrip()
        admixtures=map(str.rstrip,f.readlines())
    print newick_tree
    tree,translates=parse_newick_tree(newick_tree)
    reverse_translates={v:k for k,v in translates.items()}
    for k,c in translates.items():
        print k, ':', c
    for k,v in tree.items():
        print k,':',v
    print pretty_string(insert_children_in_tree(tree))
    print admixtures
    
def insert_admixtures(admixtures, translates, node_naming):
    admixture_proportion, _,_,_,origin, destination= admixtures.split()
    origin_key=translates[origin]
    destination_key=translates[destination]
    sink_key, source_key= node_naming(2) 
    tree[source_key]=tree[tree[origin_key][0], ]
    tree[origin_key]=[source_key, None, None,0,None]
    #tree[]=

def parse_newick_tree(newick_string):
    connects_with_lengths={}
    node_count=node()
    translates={}
    tree={}
    def reduce_newick(nstring, parent_key, node_count,tree, translates):
        nstring=nstring.strip()
        if ',' in nstring:
            i=0
            while nstring[i]!='(' and nstring[i]!=',':
                i+=1
            if nstring[i]==',':
                reduce_newick(nstring[:i], parent_key, node_count, tree,translates)
                reduce_newick(nstring[i+1:], parent_key, node_count, tree,translates)
                return None
            #i now has the index of the first (
            first_p=i
            count=0
            while count!=1:
                i+=1
                if nstring[i]=='(':
                    count-=1
                elif nstring[i]==')':
                    count+=1
            second_p=i
            first_nstring=nstring[first_p+1:second_p]
            to_come=nstring[second_p+1:]
            if to_come=='' or to_come[0]==';':
                reduce_newick(first_nstring, parent_key, node_count, tree,translates)  
            if to_come and to_come[0]==':':
                assert to_come[0]==':', 'unexptected to_come '+to_come
                node_val=node_count()
                key_name=node_val
                translates[key_name]=nstring.replace("'", '')
                reduce_newick(first_nstring, key_name, node_count,tree,translates)
                if ',' in to_come:
                    comma=to_come.split(',')
                    next_string=','.join(comma[1:])
                    reduce_newick(next_string, parent_key, node_count,tree,translates)
                    to_come=comma[0]
                blength=float(to_come.split(':')[1])
                tree[key_name]=[parent_key, None, None, blength, None]
        else:
            key,str_length= nstring.split(':')
            key_replaced=key.replace("'", '')
            tree[key_replaced]=[parent_key, None,None,float(str_length),None]
            translates[key_replaced]=nstring.replace("'", '')
        nstring.split(',')
    reduce_newick(newick_string, 'r', node_count, tree,translates)
    return tree,translates
    
if __name__=='__main__':
    read_treemix_file('sletmig/_treemix1.treeout.gz')
    