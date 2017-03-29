from Rtree_operations import get_categories, get_destination_of_lineages, propagate_married, propagate_admixtures
from copy import deepcopy
from prior import matchmake


def node_count(tree):
    '''
    Returns number of leaves, coalescence nodes, and admixture nodes in a tuple of size 3.
    '''
    
    leaves, coalescence_nodes, admixture_nodes=get_categories(tree)
    
    return len(leaves), len(coalescence_nodes), len(admixture_nodes)

def generation_counts(tree):
    '''
    returns a list of tuples of the form (no_admixtures, no_waiting_coalescences, no_sudden_coalescences, no_awaited_coalescences).
    '''
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in leaves]
    coalescences_on_hold=[]
    
    res=[]
    
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        #print 'sames', sames
        #print 'single_coalescences', single_coalescences
        #print 'admixtures', admixtures
        #print 'waiting_coalescences', waiting_coalescences
        #print 'awaited_coalescences', awaited_coalescences
        #print 'still_on_hold', still_on_hold
        res.append((len(sames), len(waiting_coalescences), len(awaited_coalescences), len(admixtures)))
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return res

def order_previously_unseparable(hierarchy, events):
    '''
    This function returns a sequence of the position of each of the events. If no ordering can be put, they are given the same number.
    '''
    pass
    

def make_dics_first_and_second(double_list):
    if double_list:
        firsts, seconds=map(list,zip(*double_list))
        dic={a:b for a,b in zip(firsts+seconds, seconds+firsts)}
        return dic, firsts, seconds
    else:
        return {},[],[]
    
    
def unique_identifier(tree):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in sorted(leaves)]
    lineages=deepcopy(ready_lineages)
    gen_to_column=range(len(ready_lineages))
    list_of_gens=[]
    coalescences_on_hold=[]
    gone=[]
    
    res=[]
    
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        res.append((len(sames), len(waiting_coalescences), len(awaited_coalescences), len(admixtures)))
        
        sames_dic, first_sames, second_sames = make_dics_first_and_second(sames)
        awaited_dic, first_awaited, second_awaited = make_dics_first_and_second(awaited_coalescences)
        waiting=waiting_coalescences.keys()
        gen=[]
        #print 'lineages',lineages
        #print 'sames', sames, sames_dic, first_sames, second_sames
        #print 'awaited', awaited_coalescences, awaited_dic, first_awaited, second_awaited
        for n,element in enumerate(lineages):
            if element in gone:
                gen.append('_')
            elif element in sames_dic:
                partner_index=lineages.index(sames_dic[element])
                if n<partner_index:
                    gen.append('c')
                else:
                    gen.append(partner_index)
            elif element in awaited_dic:
                partner_index=lineages.index(awaited_dic[element])
                if n<partner_index:
                    gen.append('c')
                else:
                    gen.append(partner_index)
            elif element in admixtures:
                gen.append('a')
            else:
                gen.append('w')
        #print 'gen',gen
        #print 'gone', gone
        list_of_gens,gone, lineages =update_lineages(list_of_gens,gen,gone, lineages, tree)
                
                
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return _list_identifier_to_string(list_of_gens)

def _list_identifier_to_string(list_of_gens):
    return '-'.join(['.'.join(map(str,[c for c in l if c!='_'])) for l in list_of_gens])

def identifier_to_tree(identifier, leaves=None, inner_nodes=None):
    pass

def update_lineages(lists, new, gone, lineages, tree):
    for n,element in enumerate(new):
        if element=='a':
            key_of_admix_child, branch_of_admix_child=lineages[n]
            admix_key=tree[key_of_admix_child][branch_of_admix_child]
            lineages.append((admix_key,1))
            lineages[n]=(admix_key,0)
        elif element=='c':
            key_of_child, branch_of_child=lineages[n]
            key=tree[key_of_child][branch_of_child]
            lineages[n]=(key,0)
        elif element=='w' or element=='_':
            pass
        else:
            gone.append(lineages[n])
    lists.append(new)
    return lists, gone, lineages
            
            
    

if __name__=='__main__':
    from Rcatalogue_of_trees import *
    
    print generation_counts(tree_good)
    print unique_identifier(tree_good)