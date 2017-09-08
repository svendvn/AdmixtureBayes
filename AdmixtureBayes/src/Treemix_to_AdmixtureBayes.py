import subprocess
from newick import parse_tree


def read_treemix_file(filename):
    reduced_filename='.'.join(filename.split(".")[:-1])
    take_copy_args=['cp', filename, filename+".tmp"]
    move_back_args=['mv', filename+'.tmp', filename]
    args=['gunzip', '-f', filename]
    subprocess.call(take_copy_args)
    subprocess.call(args)
    subprocess.call(move_back_args)
    with open(reduced_filename, 'r') as f:
        newick_tree=f.readline().rstrip()

    print type(parse_tree(newick_tree))
    t=parse_tree(newick_tree)
    print dir(t)
    print parse_tree(str(t.get_edges()[1])).get_edges()
    
if __name__=='__main__':
    read_treemix_file('sletmig/_treemix1.treeout.gz')
    