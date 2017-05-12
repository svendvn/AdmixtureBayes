from Rtree_operations import (get_categories, get_destination_of_lineages, propagate_married, 
                              propagate_admixtures, get_branch_length,update_parent_and_branch_length, 
                              get_trivial_nodes, pretty_string, insert_children_in_tree, rename_root,
                              get_admixture_proportion)
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

def get_timing(tree):
    '''
    returns a list of tuples of the form (no_admixtures, no_waiting_coalescences, no_sudden_coalescences, no_awaited_coalescences).
    '''
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in leaves]
    coalescences_on_hold=[]
    
    res={key:0.0 for key in leaves}
    count=1
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        #print 'sames', sames
        #print 'single_coalescences', single_coalescences
        #print 'admixtures', admixtures
        #print 'waiting_coalescences', waiting_coalescences
        #print 'awaited_coalescences', awaited_coalescences
        #print 'still_on_hold', still_on_hold
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))
        
        res.update({key:count for key,_ in ready_lineages})
        
        count+=1
        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            res['r']=count
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

def _list_double_to_string(list_of_doubles, digits=3):
    return '-'.join(map(str, [round(double_, digits) for double_ in list_of_doubles]))

class generate_numbered_nodes(object):
    
    def __init__(self, prefix='n', admixture_prefix='a'):
        self.count=0
        self.prefix='n'
        self.admixture_prefix='a'
        
    def __call__(self, admixture=False):
        self.count+=1
        if admixture:
            return self.admixture_prefix+str(self.count)
        return self.prefix+str(self.count)
    
class generate_predefined_list(object):
    
    def __init__(self, listi):
        self.listi=listi
        
    def __call__(self):
        return float(self.listi.pop(0))
    
def identifier_to_tree(identifier, leaves=None, inner_nodes=None, branch_lengths=None, admixture_proportions=None):
    '''
    Transforms an identifier of the form qwert-uio-asdfg-jk into a dictionary tree using the generators of leaves, inner_nodes, branch_lengths and admixture_proportions.
    '''
    
    levels=identifier.split('-')
    n_leaves=len(levels[0].split('.'))
    
    #initiate leaves
    if leaves is None:
        leaf_values=get_trivial_nodes(n_leaves)
    else:
        leaf_values=[leaves() for _ in range(n_leaves)]
    tree={leaf:[None]*5 for leaf in leaf_values}
    trace_lineages=[(leaf,0) for leaf in leaf_values]
    
    #initiate generators
    if inner_nodes is None:
        inner_nodes=generate_numbered_nodes('n')
    if branch_lengths is None:
        def f(): 
            return 1.0
        branch_lengths= f
    if admixture_proportions is None:
        def g():
            return 0.4
        admixture_proportions=g
    print branch_lengths.listi
    for level in levels:
        identifier_lineages=level.split('.')
        assert len(trace_lineages)==len(identifier_lineages), 'the number of traced lineages did not match the number of lineages in the identifier '+\
                                                               '\n\n'+'trace_lineages:'+'\n'+str(trace_lineages)+\
                                                               '\n\n'+'identifier_lineages:'+'\n'+str(identifier_lineages)
        parent_index={}
        indexes_to_be_removed=[]
        for n,identifier_lineage in enumerate(identifier_lineages):
            if identifier_lineage=='c':
                ##there is a coalecence for the n'th lineage, and it should be replaced by a new lineage
                new_key=inner_nodes()
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                tree[new_key]=[None]*5
                parent_index[n]=new_key
                trace_lineages[n]=(new_key,0)
            elif identifier_lineage=='w':
                pass
            elif identifier_lineage=='a':
                new_key=inner_nodes(admixture=True)
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                new_admixture_proportion=admixture_proportions()
                tree[new_key]=[None,None,new_admixture_proportion,None,None]
                trace_lineages[n]=(new_key,0)
                trace_lineages.append((new_key,1))
            else:
                ##there is a coalescence but this lineage disappears
                try:
                    new_key=parent_index[int(identifier_lineage)]
                except KeyError as e:
                    print e
                    print 'new_key', new_key
                    print 'parent_index', parent_index
                    print 'identifier_lineage', identifier_lineage
                    print pretty_string(insert_children_in_tree(tree))
                
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                indexes_to_be_removed.append(n)
        
        ##remove lineages
        trace_lineages=[trace_lineage for n,trace_lineage in enumerate(trace_lineages) if n not in indexes_to_be_removed]
    root_key=new_key
    del tree[root_key]
    tree=rename_root(tree, new_key)
    
    return insert_children_in_tree(tree)
              
def identifier_to_tree_clean(identifier):            
    ad2, branch_lengths, admixture_proportions=divide_triple_string(identifier)
    tree_good2= identifier_to_tree(ad2, 
                                   branch_lengths=string_to_generator(branch_lengths), 
                                   admixture_proportions=string_to_generator(admixture_proportions))
    return tree_good2

def unique_identifier_and_branch_lengths(tree):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    ready_lineages=[(key,0) for key in sorted(leaves)]
    lineages=deepcopy(ready_lineages)
    gen_to_column=range(len(ready_lineages))
    list_of_gens=[]
    coalescences_on_hold=[]
    gone=[]
    
    res=[]
    branch_lengths=[]
    admixture_proportions=[]
    
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
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
                else:
                    gen.append(partner_index)
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
            elif element in awaited_dic:
                partner_index=lineages.index(awaited_dic[element])
                if n<partner_index:
                    gen.append('c')
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
                else:
                    gen.append(partner_index)
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
            elif element in admixtures:
                gen.append('a')
                branch_lengths.append(get_branch_length(tree, element[0],element[1]))
                admixture_proportions.append(get_admixture_proportion(tree, child_key=element[0],child_branch=element[1]))
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
    return ';'.join([_list_identifier_to_string(list_of_gens),
                     _list_double_to_string(branch_lengths, 3),
                     _list_double_to_string(admixture_proportions, 3)])


def list_to_generator(listi):
    return generate_predefined_list(listi)
    
def string_to_generator(stringi, floats=True):
    if floats:
        return list_to_generator(map(str,stringi.split('-')))
    return list_to_generator(stringi.split('-'))
    
def divide_triple_string(stringbig):
    return stringbig.split(';')
    


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
    from tree_plotting import pretty_string
    
    print generation_counts(tree_good)
    print get_timing(tree_good)
    ad= unique_identifier_and_branch_lengths(tree_good)
    print pretty_string(tree_good)
    print 'identifier=',ad 
    ad2, branch_lengths, admixture_proportions=divide_triple_string(ad)
    tree_good2= identifier_to_tree(ad2, 
                                   branch_lengths=string_to_generator(branch_lengths), 
                                   admixture_proportions=string_to_generator(admixture_proportions))
    
    for key,val in tree_good2.items():
        print key,':', val
    print pretty_string(tree_good2)