from posterior import initialize_posterior
from numpy import array
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from Rtree_operations import get_pruned_tree_and_add

def get_true_posterior(wishart_file='tmp_covariance_and_multiplier.txt', 
                       true_tree_file='tmp_scaled_true_tree.txt', 
                       p=0.5, 
                       use_skewed_distr=False, 
                       wishart_df_file='tmp_wishart_DF.txt',
                       outgroup='out'):
    true_tree, nodes=read_tree_file(true_tree_file)
    print true_tree, nodes, outgroup
    x=get_pruned_tree_and_add(true_tree,outgroup)
    nodes.remove(outgroup)
    covariance, multiplier= read_wishart_file(wishart_file, nodes)
    wishart_df=read_wishart_df_file(wishart_df_file)
    posterior=initialize_posterior(covariance, M=wishart_df, p=p, use_skewed_distr=use_skewed_distr, multiplier=multiplier, nodes=nodes)
    print posterior(x)
    
def read_wishart_df_file(filename):
    with open(filename, 'r') as f:
        res=float(f.readline())
    return res    

def read_tree_file(filename):
    with open(filename, 'r') as f:
        lines=f.readlines()
        nodes=lines[0].rstrip().split()
        tree=identifier_to_tree_clean(lines[1].rstrip(), leaves=generate_predefined_list_string(nodes))
    return tree, nodes
    
def read_wishart_file(filename, nodes):
    covariance=[]
    with open(filename, 'r') as f:
        lines=f.readlines()
        nodes2=lines[0].split()
        if len(lines)==len(nodes2)+2:#in this case, there is a multiplier
            multiplier=float(lines[-1].split('=')[-1])
        else:
            multiplier=None
        for i in range(1,len(nodes2)+1):
            covariance.append(map(float,lines[i].split()[1:]))
    covariance=array(covariance)
    if nodes!=nodes2:
        mapping={val:key for key, val in enumerate(nodes2)}
        new_order=[mapping[node] for node in nodes]
    covariance=covariance[:,new_order][new_order]
    return covariance, multiplier, lines
            
if __name__=='__main__':
    print get_true_posterior()