from copy import deepcopy
from numpy import zeros, diag, ix_, outer

from tree_operations import get_number_of_admixes_on_branch
from itertools import izip

from Rtree_operations import node_is_non_admixture, node_is_coalescence, get_leaf_keys
import warnings
from __builtin__ import True

try:
    from covariance_matrix_wrapper import Covariance_Matrix2
    imported_fast_covariance=True
except ImportError:
    imported_fast_covariance=False

class Covariance_Matrix():
    
    def __init__(self, nodes_to_index):
        self.ni=nodes_to_index
        self.covmat=zeros((len(nodes_to_index), len(nodes_to_index)))
    
    def get_indices(self, nodes):
        return [self.ni[n] for n in nodes]
    
    def get_addon(self, branch_length, weights):
        return branch_length*outer(weights, weights)
    
    def update(self, branch_length, population):
        indices=self.get_indices(population.members)
        self.covmat[ix_(indices,indices)]+=self.get_addon(branch_length, population.weights)
        #self.covmat[ix_(indices,indices)]+=branch_length*outer(weights, weights)
        
    def get_matrix(self):
        return self.covmat
    
if not imported_fast_covariance:
    Covariance_Matrix2=Covariance_Matrix
    warnings.warn('Using the slow covariance matrix implemented in numpy', RuntimeWarning)
else:
    print 'Using the fast C++ implemented covariance matrix'




class Population:
    
    def __init__(self, weights,members):
        self.weights=weights
        self.members=members
        
    def get_weight(self, member):
        for m,w in zip(self.members, self.weights):
            if m==member:
                return w
        
    def subset_of_the_candidates(self,candidates):
        if any((cand in self.members for cand in candidates)):
            if any((   (cand not in self.members) or (self.get_weight(cand)<1.0) for cand in candidates)):
                return 'partly'
            else:
                return 'all'
        return 'none'
        
    def get_population_string(self, min_w, keys_to_remove=[]):
        return '.'.join(sorted([m for m,w in zip(self.members,self.weights) if w>min_w and m not in keys_to_remove]))
        
    def remove_partition(self, weight):
        #print "weight",weight
        #print "self.weight",self.weights
        n_w=[w*weight for w in self.weights]
        self.weights=[w*(1-weight) for w in self.weights]
        #print "weight",weight
        #print "self.weight",self.weights
        return Population(n_w, deepcopy(self.members))
    
    def merge_with_other(self, other):
        #print "self",self.members, self.weights
        #print "other", other.members, other.weights
     
        self.weights=[w+other.weights[other.members.index(m)] if m in other.members else w for m,w in izip(self.members,self.weights) ]
        tmpl=[(w,m) for w,m in izip(other.weights, other.members) if m not in self.members]
        if tmpl:
            x_w,x_m=zip(*tmpl)
            self.weights.extend(x_w)
            self.members.extend(x_m)
        return self
        
        #elf.pop={member:(self.pop.get(member,0.0)+other.pop.get(member,0.0)) for member in set(self.pop.keys()+other.pop.keys())}
        #print self.pop
        
    def get_contributions_as_iterable(self, branch_length):
        #print "calculating for the piece:"
        #print self.pop
        #print branch_length
        #list_of_res=[]
        for s1,w1 in self.pop.iteritems():
            for s2,w2 in self.pop.iteritems():
                #list_of_res.append(((s1,s2),w1*w2*branch_length))
                yield ((s1,s2),w1*w2*branch_length)
                #yield ((s1,s2),w1*w2*branch_length)
        #return list_of_res
        


       



    
def leave_node(key, node, population, covmat):
    if node_is_non_admixture(node): 
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch_length=node[4], population=new_pop, covmat=covmat, dependent='none')]

def follow_branch(parent_key, branch_length, population, covmat, dependent="none"):
    covmat.update(branch_length, population)
    return parent_key, population, dependent

def _add_to_waiting(dic,add,tree):
    key,pop,dep=add
    if key in dic:#this means that the node is a coalescence node
        dic[key][0][1]=pop
        dic[key][1][1]=dep
    else:
        if key=='r' or node_is_non_admixture(tree[key]):
            dic[key]=[[pop,None],[dep,"empty"]]
        else:
            dic[key]=[[pop],[dep]]
    return dic

def _full_node(key,dic):
    if key in dic:
        for dep in dic[key][1]:
            if dep=="empty":
                return False
        return True
    return False

def _merge_pops(pops):
    if len(pops)==1:
        return pops[0]
    else:
        #print pops
        return pops[0].merge_with_other(pops[1])

def _thin_out_dic(dic, taken):
    ready_nodes=[]
    #print dic
    for key,[pops, deps] in dic.items():
        #print pops, deps
        full_node=True
        for dep in deps:
            if dep is None or not (dep=="none" or _full_node(dep,dic) or dep in taken):
                full_node=False
            else:
                pass
        if full_node:
            taken.append(key)
            ready_nodes.append((key,_merge_pops(pops)))
            del dic[key]
    return dic,ready_nodes

def make_newicks(tree, node_keys=None):
    if node_keys is None:
        node_keys=sorted(get_leaf_keys(tree))
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=zip(node_keys,pops)
                
def make_covariance(tree, node_keys=None, old_cov=False):
    if node_keys is None:
        node_keys=sorted(get_leaf_keys(tree))
    #print node_keys
    #print get_leaf_keys(tree)
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=zip(node_keys,pops)
    covmat=Covariance_Matrix2({node_key:n for n,node_key in enumerate(node_keys)})
    if old_cov:
        covmat=Covariance_Matrix2({node_key:n for n,node_key in enumerate(node_keys)})
    waiting_nodes={}
    taken_nodes=[]
    while True:
        for key,pop in ready_nodes:
            upds=leave_node(key, tree[key], pop, covmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        #print 'waiting_nodes', waiting_nodes
        #print 'ready_nodes', ready_nodes
        #print 'taken_nodes', taken_nodes
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break

    return covmat.get_matrix()

class dummy_covmat(object):
    
    def update(self, *args, **kwargs):
        pass
        

def get_populations(tree, min_w=0.0, keys_to_include=None):
    
    node_keys=sorted(get_leaf_keys(tree))
    if keys_to_include is None:
        keys_to_remove=[]
    else:
        keys_to_remove=list(set(node_keys)-set(keys_to_include))
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=zip(node_keys,pops)
    waiting_nodes={}
    taken_nodes=[]
    covmat=dummy_covmat()
    pop_strings=[]
    while True:
        for key,pop in ready_nodes:
            pop_strings.append(pop.get_population_string(min_w, keys_to_remove))
            upds=leave_node(key, tree[key], pop, covmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        #print 'waiting_nodes', waiting_nodes
        #print 'ready_nodes', ready_nodes
        #print 'taken_nodes', taken_nodes
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            big_pop=ready_nodes[0][1]
            pop_strings.append(big_pop.get_population_string(min_w, keys_to_remove))
            break
    if '' in pop_strings:
        pop_strings.remove('')
    return sorted(list(set(pop_strings)))
                 
            
            
            
    


if __name__=="__main__":
    from tree_plotting import pretty_print, plot_as_directed_graph
    from Rtree_operations import insert_children_in_tree
    
    tree_clean=insert_children_in_tree({'s1':['s1s2',None, None, 0.1,None],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3':['r',None, None, 0.2, None]})
    
    tree_one_admixture=insert_children_in_tree({'s1':['s1b',None, None, 0.1,None],
          's1b':['s1s2','s3b',0.2, 0.1,0.2],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3b':['r',None, None, 0.2, None],
          's3':['s3b',None,None,0.2,None]})
    
    tree_two_admixture=insert_children_in_tree({'s1':['s1b',None, None, 0.1,None],
          's1c':['s1s2','s3b', 0.4,0.05,0.1],
          's1b':['s1c','s3a',0.2, 0.05,0.2],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3b':['r',None, None, 0.2, None],
          's3':['s3a',None,None,0.1,None],
          's3a':['s3b', None,None,0.1,None]
          })
    
    tree_two_admixture_cross=insert_children_in_tree({'s1':['s1b',None, None, 0.1,None],
          's1c':['s1s2','s3a', 0.4,0.05,0.1],
          's1b':['s1c',None,None, 0.05,None],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3b':['r','s1b', 0.4, 0.2, 0.2],
          's3':['s3a',None,None,0.1,None],
          's3a':['s3b', None,None,0.1,None]
          })
    
    tree_illegal=insert_children_in_tree({'s1':['s1b',None, None, 0.1,None],
          's1c':['s1s2','s3a', 0.4,0.05,0.1],
          's1b':['s1c','s3b',0.2, 0.05,0.2],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3b':['r',None, None, 0.2, None],
          's3':['s3a',None,None,0.1,None],
          's3a':['s3b', None,None,0.1,None]
          })
    
    tree_on_the_border=insert_children_in_tree({'s1':['c',None, None, 0.1,None],
          's2':['a',None, None,0.05,None],
          's3':['b',None,None, 0.3,None],
          'a':['b','d', 0.5,0.2,0.1],
          'c':['r','d',0.5,0.1,0.1],
          'd':['e',None,None,0.05,None],
          'b':['e',None,None,0.02,None],
          'e':['r',None,None,0.05,None]})
    
    tree_on_the_border2=insert_children_in_tree({'s1':['d',None, None, 0.1,None],
          's2':['a',None, None,0.05,None],
          's3':['e',None,None, 0.3,None],
          's4':['b',None,None, 0.3,None],
          'a':['b','c', 0.5,0.2,0.1],
          'c':['e','d',0.5,0.1,0.1],
          'b':['f',None,None,0.05,None],
          'f':['r',None,None,0.02,None],
          'e':['f',None,None,0.05,None],
          'd':['r',None,None,0.05,None]})
    
    tree_admix_to_child=insert_children_in_tree({
        's1':['r',None,None, 0.1,None],
        's2':['s2a',None,None,0.1,None],
        's3':['s3s2',None,None,0.1,None],
        's2a':['s3s2','s3s2a', 0.5,0.1,0.13],
        's3s2':['s3s2a',None,None,0.1,None],
        's3s2a':['r',None,None,0.01]
        })
    
    #print make_covariance(tree_clean,['s1','s2','s3'])
    #print make_covariance(tree_one_admixture,['s1','s2','s3'])
    print 'two admixtures same direction consistent'
    pretty_print(tree_two_admixture)
    #plot_as_directed_graph(tree_two_admixture, drawing_name='h1.BMP')
    print make_covariance(tree_two_admixture,['s1','s2','s3'])
    print 'two admixtures cross same direction'
    pretty_print(tree_illegal)
    #plot_as_directed_graph(tree_illegal, drawing_name='h2.BMP')
    print make_covariance(tree_illegal,['s1','s2','s3'])
    print 'admixture-admixture coalescent'
    #plot_as_directed_graph(tree_on_the_border, drawing_name='h3.BMP')
    pretty_print(tree_on_the_border)
    print make_covariance(tree_on_the_border,['s1','s2','s3'])
    print 'admixture on admixture'
    #plot_as_directed_graph(tree_on_the_border2, drawing_name='h4.BMP')
    pretty_print(tree_on_the_border2)
    print make_covariance(tree_on_the_border2,['s1','s2','s3','s4'])
    print 'admixture to higher in tree'
    #plot_as_directed_graph(tree_admix_to_child, drawing_name='h5.BMP')
    pretty_print(tree_admix_to_child)
    print make_covariance(tree_admix_to_child,['s1','s2','s3'])
    print 'two admixtures cross different directions'
    #plot_as_directed_graph(tree_two_admixture_cross, drawing_name='h6.BMP')
    pretty_print(tree_two_admixture_cross)
    print make_covariance(tree_two_admixture_cross,['s1','s2','s3'])
    

    from Rtree_operations import create_burled_leaved_tree
    
    from generate_prior_trees import generate_phylogeny
    from numpy import array_equal
    N=3
    def arrays_equal(r1,r2):
        same_shape= (r1.shape[0]==r2.shape[0]) and (r1.shape[1]==r2.shape[1])
        if not same_shape:
            return False
        for i in range(r1.shape[0]):
            for j in range(r2.shape[1]):
                if abs(r1[i,j]-r2[i,j]) > 1e-4:
                    print r1[i,j],r2[i,j],r1[i,j]-r2[i,j]
                    return False
        return True
    
    print get_populations(tree_on_the_border2, keys_to_include=['s1','s2','s4'])
        
    import sys
    sys.exit()
        
        
    for _ in xrange(300):
        t=generate_phylogeny(N)
        r1=make_covariance(t,['s'+str(i+1) for i in range(N)], old_cov = False)
        r2=make_covariance(t,['s'+str(i+1) for i in range(N)], old_cov = True)
        if arrays_equal(r1,r2):
            print 'arrays equal'
        else:
            print 'arrays not equal'
            print r1
            print r2
            break
        
    
        
    N=40
    tree=create_burled_leaved_tree(N,1)
    #print make_covariance(tree,['s'+str(i+1) for i in range(N)])
    def som():
        for _ in range(300):
            r=make_covariance(tree,['s'+str(i+1) for i in range(N)])
            print _
        print r
    import cProfile
     
    print cProfile.run('som()')
    
    #print cProfile.run("likewise()")
    
    
    #print make_covariance(tree_flatter_list2,["s1","s2","s3","s4"])


#         
