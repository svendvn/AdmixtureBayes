
def make_readable_leaf_name(leaf_node):
    return list(leaf_node)[0]
    
class make_read_name(object):
    
    def __init__(self):
        self.count=0
        
    def __call__(self, admix=False):
        self.count+=1
        if admix:
            return 'a'+str(self.count)
        else:
            return 'n'+str(self.count)
        

def node_structure_to_networkx(node_structure):
    edges=[]
    admixture_nodes=[]
    pure_leaves=[]
    admixture_leaves=[]
    root=[]
    coalescence_nodes=[]
    namer=make_read_name()
    string_names={}
    root_set=None
    for key,node in node_structure.items():
        plotting_type=node.get_plotting_type()
        if plotting_type=='coalescence':
            string_names[key]=namer()
            coalescence_nodes.append(string_names[key])
        elif plotting_type=='admixture':
            string_names[key]=namer(admix=True)
            admixture_nodes.append(string_names[key])
        elif plotting_type=='leaf':
            string_names[key]=make_readable_leaf_name(key)
            pure_leaves.append(string_names[key])
        elif plotting_type=='admixture_leaf':
            string_names[key]=make_readable_leaf_name(key)
            admixture_leaves.append(string_names[key])
        elif plotting_type=='root':
            string_names[key]='r'
            root.append(string_names[key])
        else:
            assert False, 'unknown plotting type'+str(plotting_type)
        
    for key, node in node_structure.items():
        parents=node.get_parents()
        for parent in parents:
            edges.append((string_names[parent.name],string_names[key]))
    return pure_leaves, admixture_leaves, coalescence_nodes, admixture_nodes, root, edges

