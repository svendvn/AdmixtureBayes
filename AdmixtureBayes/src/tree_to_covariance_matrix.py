
from copy import deepcopy
from numpy import zeros, diag, ix_, outer

from tree_operations import get_number_of_admixes_on_branch
from itertools import izip

class Population:
    
    def __init__(self, weights,members):
        self.weights=weights
        self.members=members
        
    def remove_partition(self, weight):
        n_w=[w*weight for w in self.weights]
        self.weights=[w*(1-weight) for w in self.weights]
        return Population(n_w, self.members)
    
    def merge_with_other(self, other):
        #print "merge", self.pop, other.pop
#         print "self",self.members, self.weights
#         print "other", other.members, other.weights
        self.weights=[w+other.weights[other.members.index(m)] if m in other.members else w for m,w in izip(self.members,self.weights) ]
        x_w,x_m=zip(*[(w,m) for w,m in izip(other.weights, other.members) if m not in self.members])
        self.weights.extend(x_w)
        self.members.extend(x_m)
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
    
def follow_branch(branch, population, covmat):
    branch_length=branch[1]
    if len(branch)==2:
        covmat.update(branch_length, population)
        return None, (branch[0], population), "coalesce"
    else:
        a_event=branch[2]
        a_direction=a_event[0]
        if a_direction==">":
            if len(a_event)==4:
                return branch, None, "waiting for admixture in"
            else:
                covmat.update(branch_length, population)
                population.merge_with_other(a_event[4])
                return branch[:1]+branch[3:], (None, population), "received admixture in"
        else:
            covmat.update(branch_length, population)
            a_dest=a_event[1]
            a_proportion=a_event[2]
            a_identifier=a_event[3]
            new_pop=population.remove_partition(a_proportion)
            return branch[:1]+branch[3:], ((a_dest,a_identifier),(population,new_pop)), "moving admixture out"
        
def calc_covariance_from_generators(generators, node_to_index):
    res=zeros((len(node_to_index), len(node_to_index)))
    for generator in generators:
        for (a,b), x in generator:
            res[node_to_index[a],node_to_index[b]]+=x
    return res

            
            
        
def make_covariance(tree, nodes):

    node_to_index={node:i for i,node in enumerate(nodes)}
    pops=[Population([1.0],[node]) for node in nodes]
    covmat=Covariance_Matrix(node_to_index)
    #tree=deepcopy(tree)
    tree_dict={b[1]:(b[:1]+b[2:]) for b in tree}
    tree_dict["r"]="stop"
    attaches=dict(  (branch_key,(pops[n],"ready")) for n,branch_key in enumerate(nodes) if (branch_key in nodes))

    #all_back_at_root=False
    while len(attaches)>=2:
        moves=0
        for branch_key, branch in tree_dict.iteritems():
            if branch_key in attaches and attaches[branch_key][1]=="ready" and len(attaches)>=2:
                moves+=1
                new_branch, attachments, code= follow_branch(branch, attaches[branch_key][0], covmat)
                if code == "coalesce": #the branch has coalesced with another branch
                    del attaches[branch_key]
                    new_branch_key=attachments[0]
                    new_population=attachments[1]
                    if new_branch_key in attaches:
                        new_population.merge_with_other(attaches[new_branch_key][0])
                        attaches[new_branch_key]=(new_population, "ready")
                    else:
                        attaches[new_branch_key]=(new_population, "waiting for coalescer")
                elif code == "waiting for admixture in":
                    attaches[branch_key]=(attaches[branch_key][0], "waiting for admixturer")
                elif code == "moving admixture out":
                    
                    travel_to_branch_key=attachments[0][0]
                    identifier_key=attachments[0][1]
                    new_population_stay=attachments[1][0]
                    new_population_move=attachments[1][1]
                    
                    attaches[branch_key]=(new_population_stay, "ready")
                    dest_branch=tree_dict[travel_to_branch_key]
                    for n,admix_event in enumerate(dest_branch[2::2]):
                        #print admix_event
                        if admix_event[3]==identifier_key:
                            admix_event.append(new_population_move)
                            break
                    #tree_dict[travel_to_branch_key][2+n*2]=admix_event ##CHECK IF THIS IS ACTUALLY NECESSARY
                    if travel_to_branch_key in attaches and attaches[travel_to_branch_key][1]=="waiting for admixturer" and n==0:
                        attaches[travel_to_branch_key]= (attaches[travel_to_branch_key][0],"ready")
                elif code == "received admixture in":
                    _,new_population=attachments
                    attaches[branch_key]=(new_population, "ready")
                    
                tree_dict[branch_key]=new_branch
        if moves==0:
            return -1 #this means that the tree is illegal, because it does not progress.
    return covmat.covmat




if __name__=="__main__":
    p=Population([1.0,0.5],["s2","s3"])
    #print p.pop
    #print p.remove_partition(0.5).pop    
    
    
    tree_flatter_list=[["r","s1",0.3,[">","s3",None,32123], 0.1], 
                       ["s2s3","s2",0.2],
                       ["s2s3","s3",0.15, ["<","s1",0.44,32123],0.05],
                       ["r","s2s3",0.2]]
    
    tree_flatter_list2=[["r","s1",0.3,[">","s3",None,32123], 0.1, [">","s3",None,321234],0.0], 
                       ["s2s3","s2",0.2],
                       ["s2s3","s3",0.10,["<","s1",0.1, 321234], 0.05, ["<","s1",0.44,32123],0.05],
                       ["s2s3s4", "s4", 0.9],
                       ["r","s2s3s4",0.05],
                       ["s2s3s4","s2s3",0.15]]
    
    
    
    tree_flatter_list3=[["r","s1",0.3,[">","s3",None,32123], 0.1], 
                       ["s2s3","s2",0.2],
                       ["s2s3","s3",0.15, ["<","s1",0.44,32123],0.05],
                       ["r","s2s3",0.2]]
    
    N=40
    from tree_operations import make_flat_list_no_admix
    a= make_flat_list_no_admix(N)
    
    nodes=["s"+str(i) for i in range(1,N+1)]
    print nodes
    
    print make_covariance(a,nodes)
    
    def som():
        for i in range(10):
            make_covariance(a,nodes)
        
    from numpy.random import randn
    def likewise():
        res=0
        t=randn(1000)
        for j in xrange(617370*2):
            res=t[0]*t[j%1000]
        print res
    import cProfile
     
    print cProfile.run('som()')
    
    print cProfile.run("likewise()")
    
    
    #print make_covariance(tree_flatter_list2,["s1","s2","s3","s4"])


#         