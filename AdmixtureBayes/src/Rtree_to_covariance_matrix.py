
from copy import deepcopy
from numpy import zeros, diag, ix_, outer

from tree_operations import get_number_of_admixes_on_branch
from itertools import izip

class Population:
    
    def __init__(self, weights,members):
        self.weights=weights
        self.members=members
        
    def remove_partition(self, weight):
        #print "weight",weight
        #print "self.weight",self.weights
        n_w=[w*weight for w in self.weights]
        self.weights=[w*(1-weight) for w in self.weights]
        #print "weight",weight
        #print "self.weight",self.weights
        return Population(n_w, deepcopy(self.members))
    
    def merge_with_other(self, other):
        print "self",self.members, self.weights
        print "other", other.members, other.weights
     
        self.weights=[w+other.weights[other.members.index(m)] if m in other.members else w for m,w in izip(self.members,self.weights) ]
        tmpl=[(w,m) for w,m in izip(other.weights, other.members) if m not in self.members]
        if tmpl:
            x_w,x_m=zip(*tmpl)
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


#node is p1key, p2key, adm_prop, p1length, p2length, 

def node_is_coalescence(node):
    '''
    if a node has two parents, it is not a coalescence node
    '''
    return (node[1] is None)
    
def leave_node(key, node, population, covmat):
    if node_is_coalescence(node): #if the node is coalescence it is not admixture
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat, key=key)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch_length=node[3], population=population, covmat=covmat, dependent=node[1]),
                follow_branch(parent_key=node[1],branch_length=node[4], population=new_pop, covmat=covmat, dependent=node[0])]

def follow_branch(parent_key, branch_length, population, covmat, key, dependent="none"):
    covmat.update(branch_length, population)
    return parent_key, population, dependent

def _add_to_waiting(dic,add,tree):
    key,pop,dep=add
    if key in dic:#this means that the node is a coalescence node
        dic[key][0][1]=pop
        dic[key][1][1]=dep
    else:
        if key=='r' or node_is_coalescence(tree[key]):
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
        print pops
        return pops[0].merge_with_other(pops[1])

def _thin_out_dic(dic):
    ready_nodes=[]
    print dic
    for key,[pops, deps] in dic.items():
        print pops, deps
        full_node=True
        for dep in deps:
            if dep is None or not (dep=="none" or _full_node(dep,dic)):
                full_node=False
        if full_node:
            ready_nodes.append((key,_merge_pops(pops)))
            del dic[key]
    return dic,ready_nodes
                
def make_covariance(tree, node_keys):
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=zip(node_keys,pops)
    covmat=Covariance_Matrix({node_key:n for n,node_key in enumerate(node_keys)})
    waiting_nodes={}
    while True:
        for key,pop in ready_nodes:
            upds=leave_node(key, tree[key], pop, covmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes)
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break
    return covmat
                
            
            
            
    


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
                       ["s2s3","s3",0.15, ["<","s1",0.44,32123],0.01],
                       ["r","s2s3",0.2]]
    
    tricky_tree=[['s1s2', 's1', 0.098636351184364562, 
                  ['<', 's8', 0.83683139954974017, 153099371138775232361772973428021297240L], 0.052579082611121153, 
                  ['<', 's10', 1.0285415747457356, 282660132440025697629925156001566429139L], 0.0063265189824460531, 
                  ['>', 's3', None, 27406625191307519358522646488357993197L], -0.0056101757119698756], 
                 ['s1s2', 's2', 0.040280599078741181, 
                  ['<', 's8', 0.067704508103712813, 255283802870552211589608450291393930777L], 0.036506809735954389, 
                  ['>', 's10', None, 170559928605544741257612617741920732814L], 0.12749718232452359], 
                 ['s1s2s3', 's3', 0.37919875800704428, 
                  ['<', 's1s2s3s4s5s6s7s8s9', 0.79993364364610497, 289954331610242432874138528372380781738L], 0.032213149202271731, 
                  ['<', 's1', 0.1192361491469603, 27406625191307519358522646488357993197L], 0.0036287184462564844], 
                 ['s1s2s3', 's1s2', 0.048135539422048831, 
                  ['>', 's6', None, 202660481686189436666809477261341826686L], -0.020735535843079266, 
                  ['<', 's6', 0.31457674966688715, 191177979936304796838060364226526659503L], 0.083430730505394901], 
                 ['s1s2s3s4', 's4', -0.013649273136900109, 
                  ['<', 's8', 0.61222944552619263, 22081180961291908919203915426728903768L], 0.2847373912367544], 
                 ['s1s2s3s4', 's1s2s3', 0.15600578670983092, 
                  ['>', 's1s2s3s4s5s6s7', None, 335450402234970368061483855412644502641L], 0.088439786619906008], 
                 ['s1s2s3s4s5', 's5', 0.022547719796870661, 
                  ['<', 's9', 0.38150969547757796, 59110494304151464192741970664911817463L], 0.42417976581429412],
                 ['s1s2s3s4s5', 's1s2s3s4', 0.11839880882440419], 
                 ['s1s2s3s4s5s6', 's6', 0.046019204282488416, 
                  ['>', 's1s2', None, 191177979936304796838060364226526659503L], 0.14928193215000138, 
                  ['<', 's1s2', 0.46786880853955853, 202660481686189436666809477261341826686L], 0.088202631093527648, 
                  ['>', 's1s2s3s4s5s6s7s8s9', None, 276285265485067433498037691256569293163L], 0.077264948561363025, 
                  ['<', 's7', 0.98180490305765389, 337768447685188672484787924202341384675L], 0.25815562590524782], 
                 ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.035225736930216064, 
                  ['>', 's1s2s3s4s5s6', None, 107336035275577641005667953466186666952L], 0.040231629137279651, 
                  ['>', 's1s2s3s4s5s6', None, 1168140248318658861842075078703761327L], 0.14121276979339414], 
                 ['s1s2s3s4s5s6s7', 's7', 0.10347456879334746, 
                  ['>', 's6', None, 337768447685188672484787924202341384675L], 0.23374795454940298, 
                  ['>', 's1s2s3s4s5s6s7', None, 67184406670990443088096263632659354598L], 0.075033111875982966, 
                  ['<', 's8', 0.19307617209787492, 229515252286097966861347228296660759870L], 0.23853372258498662], 
                 ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', -0.034913713180245817, 
                  ['>', 's10', None, 5629259508218995509656206860843908260L], 0.065183463466337563, 
                  ['<', 's1s2s3s4s5', 0.071193634310412909, 107336035275577641005667953466186666952L], -0.078359614659089888, 
                  ['<', 's1s2s3s4s5', 0.46729289112722089, 1168140248318658861842075078703761327L], -0.023618212058394678], 
                 ['s1s2s3s4s5s6s7s8', 's8', 0.30627958137447236, 
                  ['>', 's4', None, 22081180961291908919203915426728903768L], 0.17571755790833238, 
                  ['>', 's1', None, 153099371138775232361772973428021297240L], 0.18873329652441423, 
                  ['>', 's2', None, 255283802870552211589608450291393930777L], 0.086555822138639091, 
                  ['>', 's7', None, 229515252286097966861347228296660759870L], -0.012732814314661939], 
                 ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.093954123429733888, 
                  ['<', 's7', 0.7391847045823422, 67184406670990443088096263632659354598L], 0.032943921468016932, 
                  ['<', 's1s2s3', 0.30866771432293616, 335450402234970368061483855412644502641L], 0.046965717484200642], 
                 ['s1s2s3s4s5s6s7s8s9', 's9', 0.23977442224806264, 
                  ['>', 's5', None, 59110494304151464192741970664911817463L], 0.51019644641963857], 
                 ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', -0.041047902668444033], 
                 ['r', 's10', 0.19472363829057943, 
                  ['<', 's2', 0.12229697484604317, 170559928605544741257612617741920732814L], 0.50429995668577976, 
                  ['>', 's1', None, 282660132440025697629925156001566429139L], -0.0062763403686170961, 
                  ['<', 's1s2s3s4s5s6', 0.57530239664876304, 5629259508218995509656206860843908260L], 0.19441618198673571], 
                 ['r', 's1s2s3s4s5s6s7s8s9', 0.073232040540056165, 
                  ['>', 's3', None, 289954331610242432874138528372380781738L], 0.065984083301988475, 
                  ['<', 's6', 0.059462637145686406, 276285265485067433498037691256569293163L], 0.17723220894641167]]

    
    tricky_tree3=[['s1s2', 's1', 0.090645923291796326], 
                  ['s1s2', 's2', 0.10723905567367155], 
                  ['s1s2s3', 's3', 0.19595308363567926], 
                  ['s1s2s3', 's1s2', 0.10144196582404151], 
                  ['s1s2s3s4', 's4', 0.30915464098619166], 
                  ['s1s2s3s4', 's1s2s3', 0.085838977860683049], 
                  ['s1s2s3s4s5', 's5', 0.39764734903253096], 
                  ['s1s2s3s4s5', 's1s2s3s4', 0.099506805376701701], 
                  ['s1s2s3s4s5s6', 's6', 0.51261456406456551], 
                  ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.098414875224858062], 
                  ['s1s2s3s4s5s6s7', 's7', 0.60407587607893776], 
                  ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', 0.10442653467039886], 
                  ['s1s2s3s4s5s6s7s8', 's8', 0.70285077609461577], 
                  ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.11266789522759277], 
                  ['s1s2s3s4s5s6s7s8s9', 's9', 0.79242787177116814], 
                  ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', 0.094155170583850159], 
                  ['r', 's10', 0.89630337311324426], 
                  ['r', 's1s2s3s4s5s6s7s8s9', 0.11171881845812627]]

    tricky_tree4=[['s1s2', 's1', 0.1], 
                  ['s1s2', 's2', 0.1], 
                  ['s1s2s3', 's3', 0.2], 
                  ['s1s2s3', 's1s2', 0.1], 
                  ['s1s2s3s4', 's4', 0.30000000000000004], 
                  ['s1s2s3s4', 's1s2s3', 0.1], 
                  ['s1s2s3s4s5', 's5', 0.4], 
                  ['s1s2s3s4s5', 's1s2s3s4', 0.1], 
                  ['s1s2s3s4s5s6', 's6', 0.5], 
                  ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.1], 
                  ['s1s2s3s4s5s6s7', 's7', 0.6000000000000001], 
                  ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', 0.1], 
                  ['s1s2s3s4s5s6s7s8', 's8', 0.7000000000000001], 
                  ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.1], 
                  ['s1s2s3s4s5s6s7s8s9', 's9', 0.8], 
                  ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', 0.1], 
                  ['r', 's10', 0.9], 
                  ['r', 's1s2s3s4s5s6s7s8s9', 0.1]]
    
    tricky_tree5b=[['s1s2', 's1', 0.1], 
                   ['s1s2', 's2', 0.04039253439467916, 
                    ['<', 's7', 0.6689297039068745, 192820403950396474967875471849597410595L], 0.05960746560532085], 
                   ['s1s2s3', 's3', 0.2], 
                   ['s1s2s3', 's1s2', 0.1], 
                   ['s1s2s3s4', 's4', 0.30000000000000004], 
                   ['s1s2s3s4', 's1s2s3', 0.1], 
                   ['s1s2s3s4s5', 's5', 0.4], 
                   ['s1s2s3s4s5', 's1s2s3s4', 0.1], 
                   ['s1s2s3s4s5s6', 's6', 0.5], 
                   ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.1], 
                   ['s1s2s3s4s5s6s7', 's7', 0.4699813077288295], 
                    #['>', 's2', None, 192820403950396474967875471849597410595L, <tree_to_covariance_matrix.Population instance at 0x0000000004F42C88>], 0.13001869227117058], 
                   ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', 0.1], 
                   ['s1s2s3s4s5s6s7s8', 's8', 0.7000000000000001], 
                   ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.1], 
                   ['s1s2s3s4s5s6s7s8s9', 's9', 0.8], 
                   ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', 0.1], 
                   ['r', 's10', 0.9], 
                   ['r', 's1s2s3s4s5s6s7s8s9', 0.1]]
    
    tricky_tree5c=[['s1s2', 's1', 0.1], 
                   ['s1s2', 's2', 0.04039253439467916, 
                    ['<', 's7', 0.6689297039068745, 192820403950396474967875471849597410595L], 0.05960746560532085], 
                   ['s1s2s3', 's3', 0.2], 
                   ['s1s2s3', 's1s2', 0.1], 
                   ['s1s2s3s4', 's4', 0.30000000000000004], 
                   ['s1s2s3s4', 's1s2s3', 0.1], 
                   ['s1s2s3s4s5', 's5', 0.4], 
                   ['s1s2s3s4s5', 's1s2s3s4', 0.1], 
                   ['s1s2s3s4s5s6', 's6', 0.5], 
                   ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.1], 
                   ['s1s2s3s4s5s6s7', 's7', 0.4699813077288295, 
                    ['>', 's2', None, 192820403950396474967875471849597410595L], 0.13001869227117058], 
                   ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', 0.1], 
                   ['s1s2s3s4s5s6s7s8', 's8', 0.7000000000000001], 
                   ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.1], 
                   ['s1s2s3s4s5s6s7s8s9', 's9', 0.8], 
                   ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', 0.1], 
                   ['r', 's10', 0.9], 
                   ['r', 's1s2s3s4s5s6s7s8s9', 0.1]]


    
    tricky_tree5=[['s1s2', 's1', 0.095913830054156102], 
                  ['s1s2', 's2', 0.040673936924094066, 
                   ['<', 's7', 0.66123399429950791, 192820403950396474967875471849597410595L], 0.050567263077269393], 
                  ['s1s2s3', 's3', 0.20080279443628618], 
                  ['s1s2s3', 's1s2', 0.090604197708987061], 
                  ['s1s2s3s4', 's4', 0.31975967258982174], 
                  ['s1s2s3s4', 's1s2s3', 0.11548727176952289], 
                  ['s1s2s3s4s5', 's5', 0.40996976846681704], 
                  ['s1s2s3s4s5', 's1s2s3s4', 0.10884665546125893], 
                  ['s1s2s3s4s5s6', 's6', 0.49948421481157584], 
                  ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.0784330504240104], 
                  ['s1s2s3s4s5s6s7', 's7', 0.46193995984887365], 
                   #['>', 's2', None, 192820403950396474967875471849597410595L, <tree_to_covariance_matrix.Population instance at 0x0000000004F39F48>], 0.13205240423724723], 
                   ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', 0.086836682449862451], 
                   ['s1s2s3s4s5s6s7s8', 's8', 0.71951500182731132], 
                   ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.10397907837742325], 
                   ['s1s2s3s4s5s6s7s8s9', 's9', 0.80269345631564515], 
                   ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', 0.10465925283309137], 
                   ['r', 's10', 0.89959976967386734], 
                   ['r', 's1s2s3s4s5s6s7s8s9', 0.094523636874167108]]

    tricky_tree6=[['s1s2', 's1', 0.1], 
                  ['s1s2', 's2', 0.1], 
                  ['s1s2s3', 's3', 0.16466766863921856, 
                   ['<', 's9', 0.09453968489120956, 32402123788934016701770447745948395603L], 0.035332331360781455], 
                  ['s1s2s3', 's1s2', 0.1], 
                  ['s1s2s3s4', 's4', 0.285711961785336, 
                   ['<', 's9', 0.9875062108329644, 300501442732240316186983519279448397015L], 0.014288038214664067], 
                  ['s1s2s3s4', 's1s2s3', 0.1], 
                  ['s1s2s3s4s5', 's5', 0.4], 
                  ['s1s2s3s4s5', 's1s2s3s4', 0.1], 
                  ['s1s2s3s4s5s6', 's6', 0.5], 
                  ['s1s2s3s4s5s6', 's1s2s3s4s5', 0.1], 
                  ['s1s2s3s4s5s6s7', 's7', 0.6000000000000001], 
                  ['s1s2s3s4s5s6s7', 's1s2s3s4s5s6', 0.1], 
                  ['s1s2s3s4s5s6s7s8', 's8', 0.7000000000000001], 
                  ['s1s2s3s4s5s6s7s8', 's1s2s3s4s5s6s7', 0.1], 
                  ['s1s2s3s4s5s6s7s8s9', 's9', 0.4218560365673063, 
                   ['>', 's3', None, 32402123788934016701770447745948395603L], 0.06094517933207866, 
                   ['>', 's4', None, 300501442732240316186983519279448397015L], 0.3171987841006151], 
                  ['s1s2s3s4s5s6s7s8s9', 's1s2s3s4s5s6s7s8', 0.1], 
                  ['r', 's10', 0.9], 
                  ['r', 's1s2s3s4s5s6s7s8s9', 0.1]]


    
    N=40
    from tree_operations import make_flat_list_no_admix
    a= make_flat_list_no_admix(N)
    
    nodes=["s"+str(i) for i in range(1,N+1)]
    print nodes
    
    tree={'s1':['s1s2',None, None, 0.1,None],
          's2':['s1s2', None, None, 0.1,None],
          's1s2':['r',None, None, 0.2,None],
          's3':['r',None, None, 0.4, None]}
    
    print make_covariance(tree,['s1','s2','s3'])
    nodes2=["s"+str(i) for i in range(1,4)]
#     print make_covariance(tree_flatter_list,nodes2)
#     print make_covariance(tree_flatter_list2,nodes2)
#     print make_covariance(tree_flatter_list3,nodes2)
    
    
    

    print "before",tricky_tree6   
    #a=make_covariance(tricky_tree6,["s"+str(i) for i in range(1,11)])
    #print a
    print "after", tricky_tree6    
    #a=make_covariance(tricky_tree4,["s"+str(i) for i in range(1,11)])
    #a=make_covariance(tricky_tree,["s"+str(i) for i in range(1,11)])
    
    
    def som():
        for i in range(100):
            make_covariance(a,nodes)
        
    from numpy.random import randn
    def likewise():
        res=0
        t=randn(1000)
        for j in xrange(617370*2):
            res=t[0]*t[j%1000]
        print res
    import cProfile
     
    #print cProfile.run('som()')
    
    #print cProfile.run("likewise()")
    
    
    #print make_covariance(tree_flatter_list2,["s1","s2","s3","s4"])


#         