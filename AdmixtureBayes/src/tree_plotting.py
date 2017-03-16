from csv import writer
from Rtree_operations import to_aarhus_admixture_graph, to_networkx_format
from subprocess import call
from PIL import Image
from graphviz import Digraph
import os

file_suffix=[s+'.csv' for s in ['leaves', 'inner_nodes','edges','adm_props']]

def plot_graph(*args):
    plot_as_admixture_tree(*args)

def plot_as_admixture_tree(tree, file_prefix='', drawing_name='tmp.png', popup=True):
    aarhus_tree = to_aarhus_admixture_graph(tree)
    file_names=[file_prefix+s for s in file_suffix]
    write_aarhus_tree_to_files(aarhus_tree, file_names)
    make_R_draw_from_files(drawing_name, file_names)
    if popup:
        img=Image.open(drawing_name)
        img.show()
        
def plot_as_directed_graph(tree, file_prefix='', drawing_name='tmp.BMP', popup=True):

    leaves, admixture_nodes, coalescence_nodes, root, edges= to_networkx_format(tree)
    filename, image_format= drawing_name.split('.')
    G=Digraph('G', filename=filename)
    
    leaves_graph=Digraph('l')
    leaves_graph.node_attr.update(style='filled', color='cadetblue1')
    for leaf in leaves:
        leaves_graph.node(leaf)
        
    admixture_graph=Digraph()
    admixture_graph.node_attr.update(style='filled', fillcolor='coral1', shape='box')
    for adm_n in admixture_nodes:
        admixture_graph.node(adm_n)
        
    coalescence_graph=Digraph()
    coalescence_graph.node_attr.update(style='filled', fillcolor='greenyellow')
    for cn in coalescence_nodes:
        coalescence_graph.node(cn)
        
    G.subgraph(leaves_graph)
    G.subgraph(admixture_graph)
    G.subgraph(coalescence_graph)
    G.node('r', shape='egg', color='black', style='filled', fontcolor='white')
    G.edges(edges)
    G.format = image_format
    G.render(view=popup)
    
    
    
def write_aarhus_tree_to_files(aarhus_tree, file_names):
    for rows, name in zip(aarhus_tree, file_names):
        with open(name, "wb") as f:
            writer2 = writer(f)
            writer2.writerows(rows)
    
    
def make_R_draw_from_files(drawing_name, file_names):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cmd=['Rscript', dir_path+os.path.sep+'make_drawing.R', drawing_name]+file_names
    print cmd
    call(cmd)
    
if __name__=='__main__':
    from Rtree_operations import tree_on_the_border2, insert_children_in_tree, create_trivial_tree
    tree2=insert_children_in_tree(tree_on_the_border2)
    trouble2={'a': ['n17', 'n18', 0.5, 0.0006670327290825764, 0.04000000000000001, 's2', None], 'c': ['n15', 'r', 0.5, 0.02087163982263861, 0.4814480657456043, 'n18', None], 'n16': ['n17', None, None, 0.005272434567465561, None, 's4', 's3'], 'n17': ['n18', None, None, 0.013899593800954894, None, 'a', 'n16'], 'n15': ['r', None, None, 0.05969046586907494, None, 'c', 's1'], 's3': ['n16', None, None, 0.07815645814883887, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n15', None, None, 0.5947563021746359, None, None, None], 'n18': ['c', None, None, 0.06, None, 'a', 'n17'], 's4': ['n16', None, None, 0.00017898147838901196, None, None, None]}
    trouble3_loop={'a': ['n6', 'c', 0.5, 0.20982713110997345, 0.1, 's2', None], 'c': ['n3', 'r', 0.5, 0.00729894237428298, 0.1, 'a', None], 'f': ['n6', None, None, 0.02, None, 'n6', 's3'], 's3': ['f', None, None, 0.35, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n5', None, None, 0.09888248818230921, None, None, None], 's4': ['n3', None, None, 0.3, None, None, None], 'n3': ['n5', None, None, 0.05748492641498982, None, 'c', 's4'], 'n5': ['r', None, None, 0.036333643028417996, None, 'n3', 's1'], 'n6': ['f', None, None, 0.04017286889002655, None, 'a', 'f']}
    trouble4={'a': ['n26', 'c', 0.5, 0.027341970300003883, 0.1, 's2', None], 'c': ['n14', 'n14', 0.5, 0.6436852951194083, 0.3785128237518802, 'a', None], 'n14': ['n25', None, None, 0.6207916869626942, None, 'c', 'c'], 's3': ['r', None, None, 3.1349575990772554, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n27', None, None, 0.0, None, None, None], 's4': ['n25', None, None, 0.2760885557897766, None, None, None], 'n27': ['r', None, None, 0.0, None, 's1', 'n26'], 'n26': ['n27', None, None, 0.019998220032917055, None, 'a', 'n25'], 'n25': ['n26', None, None, 0.46155978718158974, None, 'n14', 's4']}
    
    #plot_graph({'a': ['n9', 'c', 0.5, 0.2191631972452212, 0.1, 's2', None], 'n13': ['n9', None, None, 0.09558297797101889, None, 'c', 's4'], 'c': ['n13', 'd', 0.5, 0.004417022028981122, 0.1, 'a', None], 'd': ['n10', None, None, 0.09187968643541357, None, 'c', 's1'], 'n14': ['r', None, None, 0.0, None, 's3', 'n10'], 's3': ['n14', None, None, 0.0, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['d', None, None, 0.1, None, None, None], 'n10': ['n14', None, None, 0.6792517571280886, None, 'd', 'n9'], 's4': ['n13', None, None, 1.0795285427123495, None, None, None], 'n9': ['n10', None, None, 0.9635698847039689, None, 'a', 'n13']})    
    
    #plot_graph({'a': ['n5', 'c', 0.5, 0.07182199688586655, 0.1, 's2', None], 'c': ['n2', 'r', 0.5, 0.15000000000000002, 0.15000000000000002, 'a', None], 'b': ['n4', None, None, 0.0032455783560232784, None, 'n5', 's4'], 's3': ['n5', None, None, 0.2930455087788524, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n4', None, None, 0.10695449122114763, None, None, None], 's4': ['b', None, None, 0.3, None, None, None], 'n2': ['r', None, None, 0.0596819405846814, None, 'n4', 'c'], 'n4': ['n2', None, None, 0.007072481059295323, None, 'b', 's1'], 'n5': ['b', None, None, 0.12817800311413347, None, 'a', 's3']})
    
    t3={'a': ['n37', 'c', 0.5, 1.5717637721311875, 0.1, 's1', None], 'n66': ['r', None, None, 0.008798782728668674, None, 's3', 's4'], 'c': ['n54', 'n37', 0.5, 0.771318479326775, 0.07345113788460944, 'a', None], 's3': ['n66', None, None, 0.010969920361510089, None, None, None], 's2': ['n54', None, None, 0.404441491678861, None, None, None], 's1': ['a', None, None, 0.06451508173696463, None, None, None], 's4': ['n66', None, None, 1.7305330689019498, None, None, None], 'n67': ['r', None, None, 0.24519067463109384, None, 'n54', 'n37'], 'n54': ['n67', None, None, 0.25870104556004564, None, 'c', 's2'], 'n37': ['n67', None, None, 0.9342460567572629, None, 'c', 'a']}
    t2={'a': ['n37', 'c', 0.5, 1.5717637721311875, 0.1, 's1', None], 'n66': ['r', None, None, 0.008798782728668674, None, 'n62', 's4'], 'c': ['n54', 'n37', 0.5, 0.771318479326775, 0.07345113788460944, 'a', None], 'n62': ['n66', None, None, 0.00020065386715392282, None, 'n37', 's3'], 's3': ['n62', None, None, 0.010769266494356165, None, None, None], 's2': ['n54', None, None, 0.404441491678861, None, None, None], 's1': ['a', None, None, 0.06451508173696463, None, None, None], 's4': ['n66', None, None, 1.7305330689019498, None, None, None], 'n54': ['r', None, None, 0.5038917201911395, None, 'c', 's2'], 'n37': ['n62', None, None, 0.9342460567572629, None, 'c', 'a']}
    t1={'a': ['n37', 'c', 0.5, 1.5717637721311875, 0.1, 's1', None], 'c': ['n54', 'n37', 0.5, 0.771318479326775, 0.07345113788460944, 'a', None], 'n63': ['r', None, None, 0.007697913461472917, None, 'n62', 'n54'], 'n62': ['n63', None, None, 0.008999436595822597, None, 'n37', 's3'], 's3': ['n62', None, None, 0.010769266494356165, None, None, None], 's2': ['n54', None, None, 0.404441491678861, None, None, None], 's1': ['a', None, None, 0.06451508173696463, None, None, None], 's4': ['r', None, None, 1.7305330689019498, None, None, None], 'n54': ['n63', None, None, 0.5038917201911395, None, 'c', 's2'], 'n37': ['n62', None, None, 0.9342460567572629, None, 'c', 'a']}

    #plot_graph(t1)
    #plot_graph(t2)
    #plot_graph(t3)
    
    #plot_graph({'a': ['n9992', 'r', 0.5, 0.09742527412081677, 2.3912090960430885, 's2', None], 'c': ['n9994', 'n9998', 0.5, 0.4867177654038102, 0.26803992416226596, 's3', None], 's3': ['c', None, None, 0.6407767373547782, None, None, None], 's2': ['a', None, None, 5.348010468605365, None, None, None], 's1': ['n9994', None, None, 0.06786584757783283, None, None, None], 'n10000': ['n9992', None, None, 0.02201354801127401, None, 's4', 'n9998'], 's4': ['n10000', None, None, 0.29534447309823564, None, None, None], 'n9998': ['n10000', None, None, 0.5446949533490856, None, 'c', 'n9994'], 'n9994': ['n9998', None, None, 0.29712943731271024, None, 'c', 's1'], 'n9992': ['r', None, None, 1.7634886207118632, None, 'a', 'n10000']})
    tree=create_trivial_tree(4)
    print tree
    import sys
    print sys.path
    import os
    print os.environ['PATH']
    
        
    ##IN WINDOWS YOU SHOULD PUT THE FAILED SUBPROCESS.CALL to shell=True
    plot_as_directed_graph({'n14986o': ['n14997o', None, None, 0.037596221005122325, None, 'n14994o', 'n14998n'], 'n14994o': ['n14986o', None, None, 0.5995922448996258, None, 's2', 'n14993n'], 's3': ['n14993n', None, None, 0.017991107981101026, None, None, None], 's2': ['n14994o', None, None, 0.997708485557472, None, None, None], 's1': ['n14999o', None, None, 0.5171959047036405, None, None, None], 's4': ['n14998o', None, None, 0.0011134401638378648, None, None, None], 'n14999o': ['r', None, None, 0.10861780971705044, None, 's1', 'n14999n'], 'n14999n': ['n14999o', 'r', 0.48, 0.015628043045867235, 0, 'n14997o', None], 'n14998n': ['n14986o', 'n14998o', 0.48, 0.0006314720717206017, 0, 'n14993n', None], 'n14998o': ['n14997o', None, None, 4.837436272350898e-05, None, 's4', 'n14998n'], 'n14997o': ['n14999n', None, None, 0.00379039252495712, None, 'n14986o', 'n14998o'], 'n14993n': ['n14998n', 'n14994o', 0.48, 0.00033004884190202276, 0.0, 's3', None]})
    #plot_graph(t1)
    