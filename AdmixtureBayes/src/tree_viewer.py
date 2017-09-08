from tree_plotting import plot_as_directed_graph, pretty_string
from tree_statistics import topological_identifier_to_tree_clean, identifier_to_tree_clean
from Rtree_to_coefficient_matrix import get_numbers

count=0
while True:
    var = raw_input("Please enter something: ")
    if var=='q' or var=='exit' or var=='q()' or var=='exit()':
        break
    if ';' in var:
        tree=identifier_to_tree_clean(var)
        print get_numbers(tree)
        plot_as_directed_graph(tree, drawing_name='tmp'+str(count)+'.png')
    else:
        tree=topological_identifier_to_tree_clean(var)
        print get_numbers(tree)
        plot_as_directed_graph(tree, drawing_name='tmp'+str(count)+'.png')
    count+=1