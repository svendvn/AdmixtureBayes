import numpy as np
from math import floor
from tree_statistics import identifier_to_tree, generate_predefined_list_string
from tree_plotting import plot_as_directed_graph, pretty_string
from copy import deepcopy

def read_subgraph_output(filename):
    branch_lengths=[]
    ws=[]
    with open(filename, 'r') as f:
        percentage=float(f.readline().split('=')[1])
        f.readline()#empty line
        nodes=f.readline().strip().split() #nodes line
        stree=f.readline().rstrip()
        f.readline()#empty line
        branches=f.readline().rstrip()
        while len(branches)>0:
            branch_lengths.append(map(float, branches.split()))
            branches=f.readline().rstrip()
        admixture_proportions=f.readline().rstrip()
        if len(admixture_proportions)==0:
            return percentage, np.array(branch_lengths),None
        while len(admixture_proportions)>0:
            ws.append(map(float, admixture_proportions.split()))
            admixture_proportions=f.readline().rstrip()
    return percentage, stree, nodes, np.array(branch_lengths), np.array(ws)

def get_string_confidence_interval(conf,mp):
    conf_diff=conf[1]-conf[0]
    k=np.log10(conf_diff)
    return "{0:.{prec2}f}({1:.{prec}f}-{2:.{prec}f})".format(mp, conf[0],conf[1], prec=int(-floor(k)+1), prec2=int(-floor(k)))

def make_plottable_graph(filename, conf_interval_alpha=None):
    percentage, stree, nodes, cs,ws= read_subgraph_output(filename)
    no_branch_lengths=cs.shape[1]
    if ws is None:
        all_params=cs
    else:
        all_params=np.concatenate((cs,ws),axis=1 )
    ess=calc_multiESS(all_params)
    print 'ess',ess
    marg_esss=marginal_ESSs(all_params)
    mean_params=np.mean(all_params, axis=0)
    smean_params=[str(m) for m in mean_params]
    print all_params
    if conf_interval_alpha is not None:
        x=np.percentile(all_params, q=[conf_interval_alpha/2*100, (1.0-conf_interval_alpha/2)*100], axis=0)
        print x
        print mean_params
        for i,mp in enumerate(mean_params):
            smean_params[i]=get_string_confidence_interval(x[:,i],mp)
    smean_branches=smean_params[:no_branch_lengths]
    smean_admixtures=smean_params[no_branch_lengths:]
    return (percentage,ess,(min(marg_esss),max(marg_esss))), stree, nodes, smean_branches, smean_admixtures

def plot_subgraph(filename, conf_interval_alpha=None, confidence_file=None):
    if confidence_file is None:
        confidence_filename='.'.join(filename.split('.')[:-1])+'_conf.txt'
    else:
        confidence_filename=confidence_file
    outp=make_plottable_graph(filename, conf_interval_alpha=conf_interval_alpha)
    (percentage,ess,(min_ess,max_ess)), stree, nodes, smean_branches, smean_admixtures=outp
    
    branch_generator=generate_predefined_list_string([' c'+str(i) for i in range(1,len(smean_branches)+1)])
    admix_generator=generate_predefined_list_string(['a'+str(i) for i in range(1,len(smean_admixtures)+1)])
    nodes_generator=generate_predefined_list_string(deepcopy(nodes))
    tree=identifier_to_tree(stree, leaves=nodes_generator, branch_lengths=branch_generator, admixture_proportions=admix_generator)
    print pretty_string(tree)
    plot_as_directed_graph(tree, plot_edge_lengths=True)
    with open(confidence_filename,'w') as f:
        for n,row in enumerate(smean_branches):
            f.write('c'+str(n+1)+': '+row+'\n')
        for n,row in enumerate(smean_admixtures):
            f.write('a'+str(n+1)+': '+row+'\n')    

def calc_batch_means(mat,n):
    res_mat=[]
    prev=0
    for i in range(n, mat.shape[0], n):
        res_mat.append(np.mean(mat[prev:i,],axis=0))
        prev=i
    return np.array(res_mat)

def acf(x, length=40):
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1] for i in range(1, length)])

def marginal_ESSs(mat):
    esss=[]
    for i in range(mat.shape[1]):
        rhos=acf(mat[:,i], min(50, mat.shape[0]))
        ess=1.0/(1.0+2.0*sum(rhos))*mat.shape[0]
        print 'm_ess',ess
        esss.append(ess)
    return esss
    

def calc_multiESS(mat):
    simple_cov=np.cov(mat.T)
    
    n,p=mat.shape
    bn=int(n**(0.33))
    print 'n',n
    print 'p',p
    batch_mat=calc_batch_means(mat, bn)
    print batch_mat
    batch_cov=np.cov(batch_mat.T)*bn
    print 'Lambda det=', np.linalg.det(simple_cov), 'Lambda=', simple_cov
    print 'Sigma_det=', np.linalg.det(batch_cov), 'Sigma=', batch_cov
    det_frac=np.linalg.det(simple_cov)/np.linalg.det(batch_cov)
    return det_frac**(1.0/p)*n
    
    
    
if __name__=='__main__':
    #print make_plottable_graph([0.001,0.122])
    import sys
    filname=sys.argv[1]
    plot_subgraph(filname, 0.1)
    
    