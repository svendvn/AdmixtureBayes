from Rtree_operations import (get_categories, get_destination_of_lineages, propagate_married, 
                              propagate_admixtures, get_branch_length,update_parent_and_branch_length, 
                              get_trivial_nodes, pretty_string, insert_children_in_tree, rename_root,
                              get_admixture_proportion, remove_admixture, get_admixture_proportion_from_key,
                              pretty_print, get_admixture_keys_and_proportions, tree_to_0tree,
                              direct_all_admixtures)
from copy import deepcopy
from prior import matchmake
from numpy.random import random
#import bidict
import collections
import warnings
#from tree_plotting import plot_graph, plot_as_directed_graph


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
        
        
        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            res['r']=count
            break
        
        count+=1
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
    
    
def unique_identifier(tree, leaf_order=None):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    if leaf_order is not None:
        assert set(leaves)==set(leaf_order), 'the specified leaf order did not match the leaves of the tree.'
        leaves_ordered=leaf_order
    else:
        leaves_ordered=sorted(leaves)
    ready_lineages=[(key,0) for key in leaves_ordered]
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
        for gon in gone:
            lineages.remove(gon)
        gone=[]      
                
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return _list_identifier_to_string(list_of_gens)

def admixture_sorted_unique_identifier(tree, leaf_order=None, not_opposite=True):
    '''
    because of a mis-implementation a good amount of calculations could be salvaged by the option not_opposite. It should in a bug-free environment be True.
    '''
    return unique_identifier(direct_all_admixtures(tree, not_opposite), leaf_order)

def _list_identifier_to_string(list_of_gens):
    return '-'.join(['.'.join(map(str,[c for c in l if c!='_'])) for l in list_of_gens])

def _list_double_to_string(list_of_doubles, digits=3):
    format_string='.'+str(digits)+'f'
    return '-'.join(map(str, [format(double_, format_string) for double_ in list_of_doubles]))

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
    
class generate_predefined_list_string(object):
    
    def __init__(self, listi):
        self.listi=listi
        
    def __call__(self):
        return self.listi.pop(0)
    
class same_number(object):
    
    def __init__(self, value):
        self.value=value
    
    def __call__(self):
        return self.value

class uniform_generator(object):
    
    def __call__(self):
        return random()
       
    
    
class bifacturing_tree(object):
    
    def __init__(self, string, trace_lineages):
        self.string=string
        self.trace_lineages=trace_lineages
        self.counter=0
        self.old_count=-1
        
    def pop_symbol(self):
        if self.string:
            res=self.string[0]
            if len(self.string)==1:
                remover='-'
            else:
                remover=self.string[1]
            self.string=self.string[2:]
            self.old_count= self.counter
            self.update_counter(remover)
            if self.old_count in self.trace_lineages:
                return res
            else:
                return False
    
    def update_counter(self, remover):
        if remover=='.':
            self.counter+=1
        else:
            self.counter=0
        
    def lineage_in_tracer(self):
        return (n+1 in self.trace_lineages)
    
    def update_tracer(self, key, new_value):
        self.trace_lineages[key]=new_value
        
    def remove_value_from_tracer(self, key):
        del self.trace_lineages[key]
        
    def get_lineage_count(self, key):
        return self.trace_lineages[key]
    
    def get_lineage_name(self):
        return self.trace_lineages.inv[self.old_count]
    
    def update_current_lineage_count(self, new_val):
        self.trace_lineages.forceupdate({new_val:self.old_count})
    
    def finished(self):
        return (not self.string)
    
    def get_res(self):
        assert len(self.trace_lineages)==1
        return self.trace_lineages.values()[0]
    
def non_admixture_to_newick(tree):
    leaves,_,_=get_categories(tree)
    keys_to_pops={l:l for l in leaves}
    while len(keys_to_pops)>1:
        next_gen={} #dictionary of mapping from parent to a list of children
        dup_children=[]
        dup_parents=[]
        for key in keys_to_pops:
            parent_key=tree[key][0]
            if parent_key in next_gen:
                next_gen[parent_key]=[next_gen[parent_key][0], key]
                dup_parents.append(parent_key)
                dup_children.append(next_gen[parent_key])
            else:
                next_gen[parent_key]=[key]
#         print next_gen
#         print dup_parents
#         print dup_children
        for (c1,c2),p in zip(dup_children, dup_parents):
            keys_to_pops[p]='('+','.join(sorted([keys_to_pops[c1],keys_to_pops[c2]]))+')'
            del keys_to_pops[c1]
            del keys_to_pops[c2]
    return keys_to_pops.values()[0]
    
def tree_to_prop_newicks(tree, leaves=None):
    '''
    Transforms a tree in the admixture identifier format into a dictionary of newick trees where the keys are the trees and the values are the proportion of trees of that form.
    '''      
    
    leaves,_,admixture_keys=get_categories(tree)
    
    k=len(admixture_keys)
    format_code='{0:0'+str(k)+'b}'
    
    
    n_trees={}
    max_count=65
    for i in range(2**k):
        pruned_tree = deepcopy(tree)
        bina= format_code.format(i)
        prop=1.0
        for adm_key,str_bin in zip(admixture_keys, list(bina)):
            int_bin=int(str_bin)
            if int_bin==0:
                prop*=1.0-get_admixture_proportion_from_key(tree, adm_key)
            else:
                prop*=get_admixture_proportion_from_key(tree, adm_key)
            if adm_key in pruned_tree:
                #print '------------------------------------------'
                #print 'removing', (adm_key, int_bin) , 'from tree:'
                #pretty_print(pruned_tree)
                pruned_tree=remove_admixture(pruned_tree, adm_key, int_bin)
        n_tree= non_admixture_to_newick(pruned_tree)
        if n_tree in n_trees:
            n_trees[n_tree]+=prop
        else:
            n_trees[n_tree]=prop
        if i>max_count:
            break
    return n_trees

def tree_to_mode_ntree(tree):
    '''
    Every admixture node is collapsed such that the branch with smaller weight is removed. This is not the same majority tree because
    '''
    leaves,_,admixture_keys=get_categories(tree)
    pruned_tree = deepcopy(tree)
    for adm_key in admixture_keys:
        if adm_key in pruned_tree:
            if get_admixture_proportion_from_key(tree, adm_key)>0.5:
                remove=1
            else:
                remove=0
            #print '------------------------------------------'
            #print 'removing', (adm_key, int_bin) , 'from tree:'
            #pretty_print(pruned_tree)
            pruned_tree=remove_admixture(pruned_tree, adm_key, remove)
    return non_admixture_to_newick(pruned_tree)

def tree_to_0ntree(tree):
    pruned_tree=tree_to_0tree(tree)
    return non_admixture_to_newick(pruned_tree)


def tree_to_random_ntree(tree):
    leaves,_,admixture_keys=get_categories(tree)
    pruned_tree = deepcopy(tree)
    for adm_key in admixture_keys:
        if adm_key in pruned_tree:
            #print '------------------------------------------'
            #print 'removing', (adm_key, int_bin) , 'from tree:'
            #pretty_print(pruned_tree)
            pruned_tree=remove_admixture(pruned_tree, adm_key, int(random()<0.5))
    return non_admixture_to_newick(pruned_tree)

def majority_tree(tree):
    
    all_trees= tree_to_prop_newicks(tree)
    sorted_trees=sorted(all_trees.items(), key=lambda x: x[1], reverse=True)
    return sorted_trees[0][0]
    
#      
#     levels=identifier.split('-')
#     n_leaves=len(levels[0].split('.'))
#     
#     #initiate leaves
#     if leaves is None:
#         leaf_values=get_trivial_nodes(n_leaves)
#     else:
#         leaf_values=[leaves() for _ in range(n_leaves)]
#         
#     tracer={n:leaf_val for n,leaf_val in enumerate(leaf_values)}
#     b_trees=[bifacturing_tree(identifier, tracer)]
#     
#     while any((not b_tree.finished() for b_tree in b_trees)):
#         for b_tree in b_trees:
#             if not b_tree.finished():
#                 val= b_tree.pop_symbol()
#                 if val:
#                     if val == 'c':
#                         old_lineage=b_tree.get_lineage_name()
#                         new_key=(old_lineage,None)
#                         b_tree.update_current_lineage_count(new_key)
#                     elif val=='w':
#                         pass
#                     elif val=='a':
#                         new_b_tree=deepcopy(b_tree)
#         
#         
#     for level in levels:
#         identifier_lineages=level.split('.')
#         next_states=[]
#         while states:
#             state=states.pop(0)
#             trace_lineages=state
#             assert len(trace_lineages)<=len(identifier_lineages), 'the number of traced lineages did not match the number of lineages in the identifier '+\
#                                                                    '\n\n'+'trace_lineages:'+'\n'+str(trace_lineages)+\
#                                                                    '\n\n'+'identifier_lineages:'+'\n'+str(identifier_lineages)
#             parent_index={}
#             indexes_to_be_removed=[]
#             for n,identifier_lineage in enumerate(identifier_lineages):
#                 if n in trace_lineages: # 
#                     if identifier_lineage=='c':
#                         ##there is a coalecence for the n'th lineage, and it should be replaced by a new lineage
#                         old_key=trace_lineages[n]
#                         new_key=(old_key,None)
#                         trace_lineages[n]=new_key
#                     elif identifier_lineage=='w':
#                         pass
#                     elif identifier_lineage=='a':
#                         
#                         new_key=inner_nodes(admixture=True)
#                         old_key,old_branch=trace_lineages[n]
#                         new_branch_length=branch_lengths()
#                         tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
#                         new_admixture_proportion=admixture_proportions()
#                         tree[new_key]=[None,None,new_admixture_proportion,None,None]
#                         trace_lineages[n]=(new_key,0)
#                         trace_lineages.append((new_key,1))
#                     else:
#                         ##there is a coalescence but this lineage disappears
#                         try:
#                             new_key=parent_index[int(identifier_lineage)]
#                         except KeyError as e:
#                             print e
#                             print 'new_key', new_key
#                             print 'parent_index', parent_index
#                             print 'identifier_lineage', identifier_lineage
#                             print pretty_string(insert_children_in_tree(tree))
#                         
#                         old_key,old_branch=trace_lineages[n]
#                         new_branch_length=branch_lengths()
#                         tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
#                         indexes_to_be_removed.append(n)
#             
#             ##remove lineages
#             trace_lineages=[trace_lineage for n,trace_lineage in enumerate(trace_lineages) if n not in indexes_to_be_removed]
#     root_key=new_key
#     del tree[root_key]
#     tree=rename_root(tree, new_key)
#     
#     return insert_children_in_tree(tree)
        
    
    
    
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
              
def identifier_to_tree_clean(identifier, **kwargs):            
    ad2, branch_lengths, admixture_proportions=divide_triple_string(identifier)
    tree_good2= identifier_to_tree(ad2, 
                                   branch_lengths=string_to_generator(branch_lengths), 
                                   admixture_proportions=string_to_generator(admixture_proportions),
                                   **kwargs)
    return tree_good2

def identifier_file_to_tree_clean(filename, **kwargs):
    with open(filename, 'r') as f:
        ls=f.readlines()
        if len(ls)>1 and len(ls[1])>5:
            leaves=generate_predefined_list_string(ls[0].rstrip().split())
            identifier=ls[1].rstrip()
        else:
            identifier=ls[0].rstrip()
            leaves=None
    try:
        ad2, branch_lengths, admixture_proportions=divide_triple_string(identifier)
        branch_lengths=string_to_generator(branch_lengths)
        admixture_proportions=string_to_generator(admixture_proportions)
    except ValueError as e:
        warnings.warn('Branch lengths and admixture proportions not supplied, assuming topological identifier..', UserWarning)
        ad2=identifier
        branch_lengths=uniform_generator()
        admixture_proportions=uniform_generator()
    tree_good2= identifier_to_tree(ad2, 
                                   branch_lengths=branch_lengths, 
                                   admixture_proportions=admixture_proportions,
                                   leaves=leaves, **kwargs)
    return tree_good2

def topological_identifier_to_tree_clean(identifier, **kwargs):
    tree_good2= identifier_to_tree(identifier,
                                   branch_lengths=uniform_generator(),
                                   admixture_proportions=uniform_generator(),
                                   **kwargs)
    return tree_good2

def unique_identifier_and_branch_lengths(tree, leaf_order=None):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    if leaf_order is not None:
        assert set(leaves)==set(leaf_order), 'the specified leaf order did not match the leaves of the tree.'
        leaves_ordered=leaf_order
    else:
        leaves_ordered=sorted(leaves)
    ready_lineages=[(key,0) for key in leaves_ordered]
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
                print 'entered gone!'
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
        for gon in gone:
            lineages.remove(gon)
        gone=[]
    
        #updating lineages
        coalescences_on_hold=still_on_hold+waiting_coalescences.values()
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return ';'.join([_list_identifier_to_string(list_of_gens),
                     _list_double_to_string(branch_lengths, 9),
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

def get_admixture_proportion_string(tree):
    keys, props= get_admixture_keys_and_proportions(tree)
    return '-'.join(keys)+';'+_list_double_to_string(props)



if __name__=='__main__':
    
    tree1={'s1':['n1',None,None,1.0,None,None,None],
           's2':['a1',None,None,1.0,None,None,None],
           's3':['n2',None,None,1.0,None,None,None],
           'a1':['n1','n2',0.2,0.3,0.9,'s2',None],
           'n1':['n3',None,None,1.0,None,'s1','a1'],
           'n2':['n3',None,None,1.0,None,'s3','a1'],
           'n3':['n4',None,None,0.4,None,'n2','n1'],
           's4':['n4',None,None,1.3,None,None,None],
           'n4':['r',None,None,0.3,None,'n3','s4'],
           's5':['r',None,None,1.2,None,None,None]}
    
    print unique_identifier_and_branch_lengths(tree1)
    
    from tree_warner import check
    from tree_plotting import pretty_string
    tree='c.c.w.w.a.0.1.w.w-w.a.w.c.3.w.w.w-c.w.w.0.w.w.c.9-c.w.w.w.w.0-c.w.0.w.w-c.w.0.w-c.0.w-c.0;0.012-0.038-0.05-0.062-0.029-0.068-0.068-0.097-0.005-0.048-0.083-0.024-0.123-0.057-0.024-0.047-0.058-0.144-0.079-0.04-0.035-0.032;0.086-0.219'
    tree2='c.c.w.w.a.0.1.w.w-w.a.w.c.3.w.w.w-c.w.w.0.w.w.c.6-c.w.w.w.w.0-c.w.0.w.w-c.w.0.w-c.0.w-c.0;0.012-0.038-0.05-0.062-0.029-0.068-0.068-0.097-0.005-0.048-0.083-0.024-0.123-0.057-0.024-0.047-0.058-0.144-0.079-0.04-0.035-0.032;0.586-0.219'
    tree3={'136': ['32', None, None, 0.000485164, None, 's8', 's6'], '172': ['16', None, None, 0.00101341, None, '244', '32'], '233': ['r', None, None, 0.000104723, None, 's7', None], 's9': ['32', None, None, 0.00224709, None, None, None], 's8': ['136', None, None, 0.00318685, None, None, None], 's3': ['244', None, None, 0.00230204, None, None, None], 's2': ['104', None, None, 0.00224093, None, None, None], 's1': ['103pruned', None, None, 0, None, None, None], 's7': ['233', None, None, 0.0, None, None, None], 's6': ['136', None, None, 0.0015962, None, None, None], 's5': ['52', None, None, 0.000780002, None, None, None], 's4': ['76', None, None, 0.0025597, None, None, None], '0': ['16', None, None, 0.000394219, None, 's10', '2'], '2': ['0', None, None, 0.000262294, None, 'out', '76'], '103pruned': ['104', '244', 0.875419, 0.00235895, 0, 's1', None], 'out': ['2', None, None, 0.00158848, None, None, None], 's10': ['0', None, None, 1.76418e-05, None, None, None], '244': ['172', None, None, 0.000340227, None, 's3', '103pruned'], '104': ['76', None, None, 0.000879448, None, 's2', '103pruned'], '76': ['2', None, None, 0.000144781, None, '104', 's4'], '32': ['172', None, None, 0.00148711, None, '136', 's9'], '16': ['52', None, None, 0.0041016, None, '172', '0'], '52': ['r', None, None, 0.000104723, None, '16', 's5']}
    tree4={'136': ['32', None, None, 0.000485164, None, 's8', 's6'], '172': ['16', None, None, 0.00101341, None, '244', '32'], '233': ['r', None, None, 0.000104723, None, 's7', None], 's9': ['32', None, None, 0.00224709, None, None, None], 's8': ['136', None, None, 0.00318685, None, None, None], 's3': ['244', None, None, 0.00230204, None, None, None], 's2': ['104', None, None, 0.00224093, None, None, None], 's1': ['103pruned', None, None, 0, None, None, None], 's7': ['233', None, None, 0.0, None, None, None], 's6': ['136', None, None, 0.0015962, None, None, None], 's5': ['52', None, None, 0.000780002, None, None, None], 's4': ['76', None, None, 0.0025597, None, None, None], '0': ['16', None, None, 0.000394219, None, 's10', '2'], '2': ['0', None, None, 0.000262294, None, 'out', '76'], '103pruned': ['104', '244', 0.875419, 0.00235895, 0, 's1', None], 'out': ['2', None, None, 0.00158848, None, None, None], 's10': ['0', None, None, 1.76418e-05, None, None, None], '244': ['172', None, None, 0.000340227, None, 's3', '103pruned'], '104': ['76', None, None, 0.000879448, None, 's2', '103pruned'], '76': ['2', None, None, 0.000144781, None, '104', 's4'], '32': ['172', None, None, 0.00148711, None, '136', 's9'], '16': ['52', None, None, 0.0041016, None, '172', '0'], '52': ['r', None, None, 0.000104723, None, '16', 's5']}
    ft=identifier_to_tree_clean(tree2)
    ft=tree3
    print check(ft)
    print pretty_string(ft)
    print unique_identifier_and_branch_lengths(ft)
    print tree_to_0ntree(ft)
    print tree_to_mode_ntree(ft)
    print tree_to_random_ntree(ft)
    print majority_tree(ft)
    
    from sys import exit
    exit()
    from Rcatalogue_of_trees import *
    
    from Rtree_operations import create_burled_leaved_tree
    
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
    print pretty_string(create_burled_leaved_tree(4,2))
    
    print majority_tree(tree_one_admixture)
    
