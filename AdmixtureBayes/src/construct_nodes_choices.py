from copy import deepcopy

def get_nodes(arguments, input_file, outgroup_name, reduce_node):
    
    if not arguments[0]:#this means that we should use the input file for nodes
        nodes=read_one_line(input_file)
    else:
        nodes=arguments
    before_added_outgroup=deepcopy(nodes)
    reduced_nodes=deepcopy(nodes)
    if outgroup_name in nodes:
        before_added_outgroup.remove(outgroup_name)
    else:
        nodes.append(outgroup)
    if reduce_node in reduced_nodes:      
        reduced_nodes.remove(reduce_node)
    if reduce_node not in nodes:
        nodes.append(reduce_node)
    return before_added_outgroup, nodes, reduced_nodes
        
def read_one_line(filename):
    with open(filename, 'r') as f:
        return f.readline().rstrip().split()
        