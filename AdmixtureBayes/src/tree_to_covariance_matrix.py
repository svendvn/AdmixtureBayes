
from copy import deepcopy
from numpy import zeros, diag

from tree_operations import get_number_of_admixes_on_branch

class Population:
    
    def __init__(self, weights,members):
        self.pop={m:w for m,w in zip(members,weights)}
        
    def remove_partition(self, weight):
        n_mem, n_w=zip(*[(m,w*weight) for m,w in self.pop.items()])
        self.pop={m:(w*(1-weight)) for m,w in self.pop.items()}
        return Population(n_w, n_mem)
    
    def merge_with_other(self, other):
        #print "merge", self.pop, other.pop
        self.pop={member:(self.pop.get(member,0.0)+other.pop.get(member,0.0)) for member in set(self.pop.keys()+other.pop.keys())}
        #print self.pop
        
    def get_contributions_as_iterable(self, branch_length):
        #print "calculating for the piece:"
        #print self.pop
        #print branch_length
        list_of_res=[]
        for s1,w1 in self.pop.items():
            for s2,w2 in self.pop.items():
                list_of_res.append(((s1,s2),w1*w2*branch_length))
                #yield ((s1,s2),w1*w2*branch_length)
        return list_of_res
    
def follow_branch(branch, population):
    branch_length=branch[1]
    if len(branch)==2:
        add=population.get_contributions_as_iterable(branch_length)
        return None, (branch[0], population), "coalesce", add
    else:
        a_event=branch[2]
        a_direction=a_event[0]
        if a_direction==">":
            if len(a_event)==4:
                return branch, None, "waiting for admixture in", None
            else:
                add=population.get_contributions_as_iterable(branch_length)
                population.merge_with_other(a_event[4])
                return branch[:1]+branch[3:], (None, population), "received admixture in", add
        else:
            add=population.get_contributions_as_iterable(branch_length)
            a_dest=a_event[1]
            a_proportion=a_event[2]
            a_identifier=a_event[3]
            new_pop=population.remove_partition(a_proportion)
            return branch[:1]+branch[3:], ((a_dest,a_identifier),(population,new_pop)), "moving admixture out", add
            
            
        
def make_covariance(tree, nodes):

    node_to_index={node:i for i,node in enumerate(nodes)}
    pops=[Population([1.0],[node]) for node in nodes]
    tree=deepcopy(tree)
    tree_dict={b[1]:(b[:1]+b[2:]) for b in tree}
    tree_dict["r"]="stop"
    attaches=dict(  (branch_key,(pops[n],"ready")) for n,branch_key in enumerate(nodes) if (branch_key in nodes))
    additions=[]

    #all_back_at_root=False
    while len(attaches)>=2:
        moves=0
        for branch_key, branch in tree_dict.items():
            if branch_key in attaches and attaches[branch_key][1]=="ready" and len(attaches)>=2:
                moves+=1
                new_branch, attachments, code, addition= follow_branch(branch, attaches[branch_key][0])
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
                if addition is not None:
                    additions.append(addition)
        if moves==0:
            return -1 #this means that the tree is illegal, because it does not progress.
    res=zeros((len(nodes), len(nodes)))
    for g in additions:
        for (a,b), x in g:
            #print (a,b),x
            res[node_to_index[a],node_to_index[b]]+=x
            
    return res




if __name__=="__main__":
    p=Population([1.0,0.5],["s2","s3"])
    print p.pop
    print p.remove_partition(0.5).pop    
    
    
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
    
    N=4
    from tree_operations import make_flat_list_no_admix
    a= make_flat_list_no_admix(N)
    
    nodes=["s"+str(i) for i in range(1,N+1)]
    print nodes
    
    print make_covariance(a,nodes)
    
    
    #print make_covariance(tree_flatter_list2,["s1","s2","s3","s4"])


#         