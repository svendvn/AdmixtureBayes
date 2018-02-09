import subprocess
from copy import deepcopy
#from newick import parse_tree
from Rtree_operations import pretty_string, insert_children_in_tree, insert_admixture_node_halfly, graft
from meta_proposal import new_node_naming_policy
from tree_plotting import plot_graph


class vertice_dictionary():
    
    TYPES=['Treemix_V',
           'Treemix_N',
           'AdmB']
    
    def __init__(self):
        
        self.treemix_numbered_vertices=[]
        self.treemix_newick_format=[]
        self.admixture_string_nodes=[]
        
        
        
    def get_list(self, type):
        if type==vertice_dictionary.TYPES[0]:
            return self.treemix_numbered_vertices
        elif type==vertice_dictionary.TYPES[1]:
            return self.treemix_newick_format
        elif type==vertice_dictionary.TYPES[2]:
            return self.admixture_string_nodes
        else:
            assert False, 'wrong type provided'
    
    def get_3_lists(self, type1, type2, type3=None):
        available_types=deepcopy(vertice_dictionary.TYPES)
        list1=self.get_list(type1)
        available_types.remove(type1)
        list2=self.get_list(type2)
        available_types.remove(type2)
        list3=self.get_list(available_types[0])
        return list1,list2,list3
            
            
    def insert_mapping(self, from_object, to_object, from_type, to_type):
        from_list, to_list, third_list=self.get_3_lists(from_type, to_type)
        index=None
        if from_object in from_list:
            index=from_list.index(from_object)
        if to_object in to_list:
            index=to_list.index(to_object)
        if index is None:
            from_list.append(from_object)
            to_list.append(to_object)
            third_list.append(None)
        else:
            #print 'index',index
            if from_object in from_list:
                to_list[index]=to_object
            if to_object in to_list:
                from_list[index]=from_object

    def get_value(self, from_object, from_type, to_type):
        from_list, to_list, _ = self.get_3_lists(from_type, to_type)
        index=from_list.index(from_object)
        return to_list[index]
    
    
    
    def __str__(self):
        res=''
        if not self.treemix_newick_format:
            return res
        for adm, treemix_V, treemix_N in zip(self.admixture_string_nodes, self.treemix_numbered_vertices, self.treemix_newick_format):
            res+='{:6}'.format(str(adm))
            res+='{:6}'.format(str(treemix_V))
            res+=str(treemix_N)+'\n'
        return res[:-1] #erasing the last line change
        
        
class node(object):
    
    def __init__(self):
        self.count=0
        
    def __call__(self):
        self.count+=1
        print 'n'+str(self.count)
        return 'n'+str(self.count)
    
class adm_node(object):
    
    def __init__(self):
        self.count=0
        
    def __call__(self):
        self.count+=1
        return 'a'+str(self.count), 'x'+str(self.count)
    

def unzip_file(filename):
    reduced_filename='.'.join(filename.split(".")[:-1])
    take_copy_args=['cp', filename, filename+".tmp"]
    move_back_args=['mv', filename+'.tmp', filename]
    args=['gunzip', '-f', filename]
    subprocess.call(take_copy_args)
    subprocess.call(args)
    subprocess.call(move_back_args)
    return reduced_filename

def match_vertices(filename_vertices, vd):
    admixture_vertices=[]
    with open(filename_vertices, 'r') as f:
        for lin in f.readlines():
            a=lin.split()
            is_root=(a[2]=='ROOT')
            is_mig=(a[3]=='MIG')
            if not is_root and not is_mig:
                V_name=a[0]
                N_name=a[-1].rstrip()
                vd.insert_mapping(from_object=N_name, to_object=V_name, from_type='Treemix_N', to_type='Treemix_V')
            if is_mig:
                admixture_vertices.append((a[0],(a[5],a[6])))
    return vd, admixture_vertices
    
def get_edge_lengths(filename_edges):
    edges={}
    with open(filename_edges, 'r') as f:
        for lin in f.readlines():
            a=lin.split()
            edges[(a[0],a[1])]=a[2]
    return edges

def parse_admixtures(admixture_strings):
    from_to_weights={}
    for adm_string in admixture_strings:
        a=adm_string.split()
        weight=float(a[0])
        from_N=a[4]
        to_N=a[5]
        from_to_weights[from_N]=(to_N,weight)
    return from_to_weights

        
def add_admixtures(tree, vd, adm_vertices, edges, admixtures):
    a_names=adm_node()
    for source_key_V, (source_parent_key_V, source_child_key_V) in adm_vertices:
        source_child_key_B=vd.get_value(source_child_key_V, 'Treemix_V','AdmB')
        source_child_key_N=vd.get_value(source_child_key_V, 'Treemix_V','Treemix_N')
        sink_child_key_N, weight=admixtures[source_child_key_N]
        sink_child_key_B=vd.get_value(sink_child_key_N, 'Treemix_N', 'AdmB')
        sink_child_V=vd.get_value(sink_child_key_N, 'Treemix_N', 'Treemix_V')
        sink_name, source_name=a_names()
        
        t3=float(edges[(source_key_V, source_child_key_V)])
        t4=float(edges[(source_parent_key_V, source_key_V)])
        u1=t3/(t3+t4)
        
        
        tree=insert_admixture_node_halfly(tree, sink_child_key_B, 0, insertion_spot=0.000001, admix_b_length=0.000001, new_node_name=sink_name, admixture_proportion= 1-weight)
        #print 'tree after inserting admixture', tree
        tree=graft(tree, sink_name, source_child_key_B, u1, source_name, 0, remove_branch=1)
    return tree


def read_treemix_file(filename_treeout, filename_vertices, filename_edges):
    np=new_node_naming_policy()
    if filename_treeout.endswith('.gz'):
        filename_treeout=unzip_file(filename_treeout)
    if filename_vertices.endswith('.gz'):
        filename_vertices=unzip_file(filename_vertices)
    if filename_edges.endswith('.gz'):
        filename_edges=unzip_file(filename_edges)
    with open(filename_treeout, 'r') as f:
        newick_tree=f.readline().rstrip()
        admixtures=parse_admixtures(map(str.rstrip,f.readlines()))
        
    #print newick_tree
    tree,translates=parse_newick_tree(newick_tree)
    vd=vertice_dictionary()
    for adm_key, treemix_N_key in translates.items():
        vd.insert_mapping(adm_key, treemix_N_key, 'AdmB', 'Treemix_N')
    #print '-------------------------'
    #print vd
    vd, adm_vertices=match_vertices(filename_vertices, vd)
    #matched_admixtures=match_admixtures(admixtures, adm_vertices)
   # print '-------------------------'
   # print vd
   # print adm_vertices
    edges=get_edge_lengths(filename_edges)
  #  print edges
    tree=insert_children_in_tree(tree)
    reverse_translates={v:k for k,v in translates.items()}
#     for k,c in translates.items():
#         print k, ':', c
#     for k,v in tree.items():
#         print k,':',v
#     print translates
#     print admixtures
    tree=add_admixtures(tree, vd, adm_vertices, edges, admixtures)
    return tree
    
    
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
                translates[key_name]=("("+first_nstring+")"+to_come.split(',')[0]).replace("'", '')
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
        #nstring.split(',')
    reduce_newick(newick_string, 'r', node_count, tree,translates)
    return tree,translates
    
if __name__=='__main__':
    tree=read_treemix_file('../../../../Dropbox/Bioinformatik/AdmixtureBayes/treemix_example3/new_one2.treeout.gz',
                           '../../../../Dropbox/Bioinformatik/AdmixtureBayes/treemix_example3/new_one2.vertices',
                           '../../../../Dropbox/Bioinformatik/AdmixtureBayes/treemix_example3/new_one2.edges')
    #plot_graph(tree)
    import numpy as np
    print pretty_string(tree)
    from Rtree_to_covariance_matrix import make_covariance
    from reduce_covariance import reduce_covariance, Areduce
    cov=make_covariance(tree, node_keys=['out']+['s'+str(i) for i in range(1,10)])
    cov2=np.loadtxt( '../../../../Dropbox/Bioinformatik/AdmixtureBayes/treemix_example3/anew.txt')
    np.set_printoptions(precision=6,  linewidth=200, suppress=True)
    print cov-cov2
    print reduce_covariance(cov-cov2,0)
    print Areduce(cov-cov2)
    
    