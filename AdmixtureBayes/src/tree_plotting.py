from csv import writer
from Rtree_operations import to_aarhus_admixture_graph
from subprocess import call
from PIL import Image
import os

file_suffix=[s+'.csv' for s in ['leaves', 'inner_nodes','edges','adm_props']]

def plot_graph(tree, file_prefix='', drawing_name='tmp.png', popup=True):
    aarhus_tree = to_aarhus_admixture_graph(tree)
    file_names=[file_prefix+s for s in file_suffix]
    write_aarhus_tree_to_files(aarhus_tree, file_names)
    make_R_draw_from_files(drawing_name, file_names)
    if popup:
        img=Image.open(drawing_name)
        img.show()
    
def write_aarhus_tree_to_files(aarhus_tree, file_names):
    for object, name in zip(aarhus_tree, file_names):
        with open(name, "wb") as f:
            writer2 = writer(f)
            writer2.writerows(object)
    
    
def make_R_draw_from_files(drawing_name, file_names):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cmd=['Rscript', dir_path+os.path.sep+'make_drawing.R', drawing_name]+file_names
    print cmd
    call(cmd, shell=True)
    
if __name__=='__main__':
    from Rtree_operations import tree_on_the_border2, insert_children_in_tree
    tree2=insert_children_in_tree(tree_on_the_border2)
    plot_graph(tree2)
    