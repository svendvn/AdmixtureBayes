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
    call(cmd)
    
if __name__=='__main__':
    from Rtree_operations import tree_on_the_border2, insert_children_in_tree
    tree2=insert_children_in_tree(tree_on_the_border2)
    trouble2={'a': ['n17', 'n18', 0.5, 0.0006670327290825764, 0.04000000000000001, 's2', None], 'c': ['n15', 'r', 0.5, 0.02087163982263861, 0.4814480657456043, 'n18', None], 'n16': ['n17', None, None, 0.005272434567465561, None, 's4', 's3'], 'n17': ['n18', None, None, 0.013899593800954894, None, 'a', 'n16'], 'n15': ['r', None, None, 0.05969046586907494, None, 'c', 's1'], 's3': ['n16', None, None, 0.07815645814883887, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n15', None, None, 0.5947563021746359, None, None, None], 'n18': ['c', None, None, 0.06, None, 'a', 'n17'], 's4': ['n16', None, None, 0.00017898147838901196, None, None, None]}
    trouble3_loop={'a': ['n6', 'c', 0.5, 0.20982713110997345, 0.1, 's2', None], 'c': ['n3', 'r', 0.5, 0.00729894237428298, 0.1, 'a', None], 'f': ['n6', None, None, 0.02, None, 'n6', 's3'], 's3': ['f', None, None, 0.35, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n5', None, None, 0.09888248818230921, None, None, None], 's4': ['n3', None, None, 0.3, None, None, None], 'n3': ['n5', None, None, 0.05748492641498982, None, 'c', 's4'], 'n5': ['r', None, None, 0.036333643028417996, None, 'n3', 's1'], 'n6': ['f', None, None, 0.04017286889002655, None, 'a', 'f']}
    trouble4={'a': ['n26', 'c', 0.5, 0.027341970300003883, 0.1, 's2', None], 'c': ['n14', 'n14', 0.5, 0.6436852951194083, 0.3785128237518802, 'a', None], 'n14': ['n25', None, None, 0.6207916869626942, None, 'c', 'c'], 's3': ['r', None, None, 3.1349575990772554, None, None, None], 's2': ['a', None, None, 0.05, None, None, None], 's1': ['n27', None, None, 0.0, None, None, None], 's4': ['n25', None, None, 0.2760885557897766, None, None, None], 'n27': ['r', None, None, 0.0, None, 's1', 'n26'], 'n26': ['n27', None, None, 0.019998220032917055, None, 'a', 'n25'], 'n25': ['n26', None, None, 0.46155978718158974, None, 'n14', 's4']}
    plot_graph({'a': ['n30', 'c', 0.5, 0.2761275660998635, 0.1, 's4', None], 'n48': ['n47', None, None, 0.0, None, 's1', 's2'], 'c': ['n30', 'n43', 0.5, 0.01948983019409606, 2.983732832048786, 'a', None], 's3': ['n43', None, None, 0.018889505288668066, None, None, None], 's2': ['n48', None, None, 0.0016763695708895115, None, None, None], 'n43': ['r', None, None, 1.4916512302504947, None, 'c', 's3'], 'n47': ['r', None, None, 0.0, None, 'n48', 'n30'], 's4': ['a', None, None, 0.12666162681264667, None, None, None], 'n30': ['n47', None, None, 0.6554252923228383, None, 'c', 'a'], 's1': ['n48', None, None, 0.0, None, None, None]})
    