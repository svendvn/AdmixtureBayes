from tree_plotting import plot_as_directed_graph, plot_graph
from Rproposal_admix import addadmix, deladmix
from Rproposal_regraft import make_regraft
from Rproposal_rescale import rescale
from Rtree_operations import create_trivial_tree, make_consistency_checks, get_number_of_admixes, get_trivial_nodes
from numpy.random import choice
from time import sleep as wait

def _get_new_nodes(i,k):
    if k==2:
        return ['x'+str(i)+a for a in ['a','b']]
    else:
        return 'x'+str(i)

def topological_support(start_tree, n=10000, nodes=None):
    tree=start_tree
    score=0
    for i in range(n):
        prop_index=choice(3,1)
        if prop_index==0 and get_number_of_admixes(tree)>0:
            new_tree=deladmix(tree)
            score-=1
        elif prop_index==1:
            new_tree=make_regraft(tree, _get_new_nodes(i, prop_index))
        elif prop_index==2:
            new_tree=addadmix(tree, _get_new_nodes(i, prop_index))
            score+=1
        else:
            new_tree=tree
        consistent, information = make_consistency_checks(new_tree,nodes)
        if consistent:
            if get_number_of_admixes(new_tree)<50:
                tree=new_tree
        else:
            print information
            print 'last_good_tree', tree
            print 'new_bad_tree', new_tree
            plot_as_directed_graph(tree, drawing_name='before.png')
            wait(1)
            plot_as_directed_graph(new_tree, drawing_name='after.png')
            wait(1)
        if i%1000==0:
            plot_as_directed_graph(tree, drawing_name='number_'+str(i)+'K.png')
            wait(1)
    plot_as_directed_graph(tree, drawing_name='final_tree.png')
    wait(1)
    return score

def proposal_support(start_tree, n=10000, nodes=None):
    tree=start_tree
    score=0
    for i in range(n):
        prop_index=choice(4,1)
        if prop_index==0 and get_number_of_admixes(tree)>0:
            new_tree=deladmix(tree)
            score-=1
        elif prop_index==1:
            new_tree=make_regraft(tree, _get_new_nodes(i, prop_index))
        elif prop_index==2:
            new_tree=addadmix(tree, _get_new_nodes(i, prop_index))
            score+=1
        elif prop_index==3:
            new_tree=rescale(tree)
        else:
            new_tree=tree
        consistent, information = make_consistency_checks(new_tree,nodes)
        if consistent:
            if get_number_of_admixes(new_tree)<50:
                tree=new_tree
        else:
            print information
            print 'last_good_tree', tree
            print 'new_bad_tree', new_tree
            plot_as_directed_graph(tree, drawing_name='before.png')
            wait(1)
            plot_as_directed_graph(new_tree, drawing_name='after.png')
            wait(1)
            break
        if i%1000==0:
            print tree
            plot_as_directed_graph(tree, drawing_name='number_'+str(i)+'K.png')
            wait(1)
    plot_as_directed_graph(tree, drawing_name='final_tree.png')
    wait(1)
    return score

if __name__=='__main__':
    s_tree=create_trivial_tree(15)
    plot_as_directed_graph(s_tree)
    wait(1)
    print proposal_support(s_tree, nodes=get_trivial_nodes(15))
    