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
        print 'sames', sames
        print 'single_coalescences', single_coalescences
        print 'admixtures', admixtures
        print 'waiting_coalescences', waiting_coalescences
        print 'awaited_coalescences', awaited_coalescences
        print 'still_on_hold', still_on_hold
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
    
    
    
    
def unique_identifier(tree):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in sorted(leaves)]
    live_lineages=deepcopy(ready_lineages)
    gen_to_column=range(len(ready_lineages))
    list_of_gens=[]
    coalescences_on_hold=[]
    gone=[]
    
    res=[]
    
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        res.append((len(sames), len(waiting_coalescences), len(awaited_coalescences), len(admixtures)))
        
        sames_dic={a:b for a,b in sames}.update({b:a for a,b in sames})
        first_sames, second_sames=list(*sames)
        awaited_dic={a:b for a,b in awaited_coalescences}.update({b:a for a,b in awaited_coalescences})
        first_awaited, second_awaited=list(*awaited_coalescences)
        waiting=waiting_coalescences.keys()
        gen=[]
        for n,element in live_lineages:
            if element in gone:
                gen.append('_')
            if element in first_sames:
                gen.append('c')
            elif element in second_sames:
                gen.append(live_lineages.index(sames_dic[element]))
            elif element in first_awaited:
                gen.append('c')
            elif element in second_awaited:
                gen.append(live_lineages.index(awaited_dic[element]))
            elif element in admixtures:
                gen.append('a')
            else:
                gen.append('w')
        list_of_gens, gen_to_column =update_lineages(list_of_gens,gen, gen_to_column)
                
                
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return res

def update_lineages(lists, new, gen_to_column):
    new_gen_to_column
    for n,(column, element) in zip(new, gen_to_column):
        if element=='a':
            n
            
            
    

if __name__=='__main__':
    from Rcatalogue_of_trees import *
    
    print generation_counts(tree_good)