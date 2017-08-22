from tree_plotting import plot_as_directed_graph, pretty_string
from tree_statistics import topological_identifier_to_tree_clean, identifier_to_tree_clean

count=0
while True:
    var = raw_input("Please enter something: ")
    if var=='q' or var=='exit' or var=='q()' or var=='exit()':
        break
    if ';' in var:
        print pretty_string(identifier_to_tree_clean(var))
        plot_as_directed_graph(identifier_to_tree_clean(var), drawing_name='tmp'+str(count)+'.png')
    else:
        plot_as_directed_graph(topological_identifier_to_tree_clean(var), drawing_name='tmp'+str(count)+'.png')
    count+=1