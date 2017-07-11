from tree_statistics import identifier_to_tree_clean



def get_covariance(stages_to_go_through=[], input, nodes=None):
    

def read_tree(input, nodes):
    if isinstance(input, basestring):
            if not ';' in input:
                input=read_one_line(filename=input)
            return identifier_to_tree_clean(input, nodes=nodes)
        else:
            return input

def read_input(stage, input, nodes):
    if stage==1:
        return input
    if stage==2:
        return input
    if stage==3:
        return read_tree(input, nodes)
    if stage==4:
        return read_tree(input, nodes)
    if stage==5:
        return input
    if stage==6:
        
def read_covariance_matrix(filename):
    
    
        
def read_one_line(filename):
    with open(filename, 'r') as f:
        return f.readline().rstrip()