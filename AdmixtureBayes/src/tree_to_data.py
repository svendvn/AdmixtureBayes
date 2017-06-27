from Rtree_to_covariance_matrix import make_covariance
from scipy.stats import wishart
from Rtree_operations import find_rooted_nodes, get_real_parents, pretty_string, get_no_leaves, node_is_non_admixture, node_is_admixture, node_is_leaf_node, node_is_coalescence, get_real_children_root, get_trivial_nodes
from tree_statistics import get_timing
import subprocess
from numpy import loadtxt, cov, array, mean, vstack, sum, identity, insert, hstack, vsplit, amin
from numpy.linalg import det
from copy import deepcopy
from load_data import read_data


def tree_to_data_perfect_model(tree, df):
    m=make_covariance(tree)
    r=m.shape[0]
    m=wishart.rvs(df=r*df-1, scale=m/(r*df))
    return m

def normalise(m):
    return m-max(0,amin(m))

def file_to_emp_cov(filename, reduce_column=None):
    
    dat=[]
    with open(filename, 'r') as f:
        f.readline()
        for l in f.readlines():
            print l
            dat.append(map(float, l.split()[1:]))
    m=array(dat)
    if reduce_column is not None:
        m=reduce_covariance(m, reduce_column)
        m=normalise(m)
    return m

def supplementary_text_ms_string():
    return 'ms 400 400 -t 200 -r 200 500000 -I 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 \
20 20 20 20 20 -en 0.00270 20 0.025 -ej 0.00275 20 19 -en 0.00545 19 0.025 -ej 0.00550 \
19 18 -en 0.00820 18 0.025 -ej 0.00825 18 17 -en 0.01095 17 0.025 -ej 0.011 17 16 -en \
0.01370 16 0.025 -ej 0.01375 16 15 -en 0.01645 15 0.025 -ej 0.01650 15 14 -en 0.01920 \
14 0.025 -ej 0.01925 14 13 -en 0.02195 13 0.025 -ej 0.02200 13 12 -en 0.02470 12 0.025 \
-ej 0.02475 12 11 -en 0.02745 11 0.025 -ej 0.02750 11 10 -en 0.03020 10 0.025 -ej 0.03025 \
10 9 -en 0.03295 9 0.025 -ej 0.03300 9 8 -en 0.03570 8 0.025 -ej 0.03575 8 7 -en 0.03845 \
7 0.025 -ej 0.03850 7 6 -en 0.04120 6 0.025 -ej 0.04125 6 5 -en 0.04395 5 0.025 -ej \
0.04400 5 4 -en 0.04670 4 0.025 -ej 0.04675 4 3 -en 0.04945 3 0.025 -ej 0.04950 3 2 \
-en 0.05220 2 0.025 -ej 0.05225 2 1'


def tree_to_ms_command(rtree, sample_per_pop=50, nreps=2, 
                       theta=0.4, sites=500000, recomb_rate=1):
    tree=deepcopy(rtree)
    if recomb_rate is None:
        rec_part=' -s '+str(sites)
    else:
        rec_part=' -r '+str(recomb_rate)+ ' '+str(sites)
    n=get_no_leaves(tree)
    callstring='ms '+str(sample_per_pop*n)+' '+str(nreps)+' -t '+ str(theta)+' ' +rec_part + ' '
    callstring+=' -I '+str(n)+' '+' '.join([str(sample_per_pop) for _ in xrange(n)])+' '
    times=get_timing(tree)
    max_v=max(times.values())
    times={k:v for k,v in times.items()}
    print times
    tree=extend_branch_lengths(tree,times)
    print pretty_string(tree)
    callstring+=construct_ej_en_es_string(tree, times, n)
    return callstring
    
def extend_branch_lengths(tree, times):
    for key, node in tree.items():
        for n,parent_key in enumerate(get_real_parents(node)):
            pseudo_time=times[parent_key]-times[key]
            tree[key][3+n]=(tree[key][3+n], pseudo_time)
    return tree



def construct_ej_en_es_string(tree, times, no_leaves):
    s_times=sorted(zip(times.values(),times.keys()))
    dic_of_lineages={(key,0):(n+1) for n,(time,key) in enumerate(s_times[:no_leaves])}
    population_count=len(dic_of_lineages)
    res_string=''
    for time,key in s_times:
        if key=='r':
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
            break
        node=tree[key]
        if node_is_leaf_node(node):
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_coalescence(node):
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=calculate_pop_size(node[3])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size)+' '
        if node_is_admixture(node):
            population_count+=1
            i=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))[0]
            res_string+='-es '+str(time)+' '+str(i)+' '+str(1.0-node[2])+' '
            dic_of_lineages[(key,0)]=i
            dic_of_lineages[(key,1)]=population_count
            pop_size1=calculate_pop_size(node[3])
            pop_size2=calculate_pop_size(node[4])
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(pop_size1)+' '
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,1)])+' '+str(pop_size2)+' '
    return res_string
            
        
def get_affected_populations(dic_of_lineages, children_branches):
    return [dic_of_lineages[children_branch] for children_branch in children_branches] 
    
def calculate_pop_size(tup):
    drift, actual=tup
    return actual/drift

def call_ms_string(ms_string, sequence_file):
    with open(sequence_file, 'w') as f:
        p = subprocess.Popen(ms_string, stdout=subprocess.PIPE, shell=True)
        line_number = 0
        for line in p.stdout.readlines():
            line_number += 1
            if line_number >= 5 and line and (line[0]=='0' or line[0]=='1'):
                #print len(line)
                print line_number, line[:4]
                f.write(line)
            else:
                print line_number,':', line.rstrip()
    return 0

def calculate_covariance_matrix(file='tmp.txt', samples_per_pop=20, no_pops=4, n_reps=1):
    data=[]
    with open(file, 'r') as f:
        for r in f.readlines():
            print r[:5]
            data.append(map(int,list(r.rstrip())))
    m= array(data)
    if n_reps>1:#reorder the data so that there are more SNPs in stead of more samples/populations
        m=hstack(vsplit(m, n_reps))
    #print 'm',m
    total_mean=mean(m, axis=0)
    #print 'total_mean', total_mean
    ps=tuple([mean(m[(i*samples_per_pop):((i+1)*samples_per_pop), : ], axis=0)-total_mean for i in xrange(no_pops)])
    p=vstack(ps)
    #print 'p',p
    #print p.shape
    #print p.dot(p.T)/(p.shape[1])
    return p.dot(p.T)/(p.shape[1])

def reduce_covariance(covmat, subtracted_population_index):
    reducer=insert(identity(covmat.shape[0]-1), subtracted_population_index, -1, axis=1)
    return reducer.dot(covmat).dot(reducer.T)

def ms_to_treemix2(filename='tmp.txt', samples_per_pop=20, no_pops=4, n_reps=1, filename2='tmp.treemix_in'):
    with open(filename, 'r') as f:
        with open(filename2, 'w') as e:
            e.write(' '.join(get_trivial_nodes(no_pops))+'\n')
            pop_count=0
            rep_count=0
            count=0
            data=[]
            s_vecs=[]
            for r in f.readlines():
                
                data.append(map(int,list(r.rstrip())))
                count+=1
                
                if count==samples_per_pop:
                    
                    count=0
                    pop_count+=1
                    
                    s_vec=sum(array(data), axis=0)
                    s_vecs.append(s_vec)
                    
                    data=[]
                    
                    print rep_count, pop_count
                    
                    if pop_count==no_pops:
                        
                        pop_count=0
                        rep_count+=1
                        
                        for s in zip(*s_vecs):
                            e.write(' '.join([str(a)+','+str(samples_per_pop-a) for a in s])+'\n')
                        
                        s_vecs=[]
                        
                        if rep_count>=n_reps:
                            break

    filename2_gz=filename2+'.gz'
    subprocess.call(['gzip','-f','-k', filename2])
    return read_data(filename2_gz, blocksize=1000 ,outgroup='s1', noss=False, nodes=get_trivial_nodes(no_pops))
    
def ms_to_treemix(filename='tmp.txt', samples_per_pop=20, no_pops=4, n_reps=1, filename2='tmp.treemix_in'):
    data=[]
    with open(filename, 'r') as f:
        for r in f.readlines():
            print r[:5]
            data.append(map(int,list(r.rstrip())))
    m= array(data)
    if n_reps>1:#reorder the data so that there are more SNPs in stead of more samples/populations
        print m.shape
        print 'samples, pops, reps', samples_per_pop, no_pops, n_reps
        m=hstack(vsplit(m, n_reps))
        
    #print samples_per_pop
    sums=tuple([sum(m[(i*samples_per_pop):((i+1)*samples_per_pop), : ], axis=0) for i in xrange(no_pops)])
    #print sums, 'sums'
    with open(filename2, 'w') as f:
        f.write(' '.join(get_trivial_nodes(no_pops))+'\n')
        for s_vec in zip(*sums):
            f.write(' '.join([str(s)+','+str(samples_per_pop-s) for s in s_vec])+'\n')
    filename2_gz=filename2+'.gz'
    subprocess.call(['gzip','-f','-k', filename2])
    return read_data(filename2_gz, blocksize=10000 ,outgroup='s3', noss=True)

def trace_from_root(tree, init_freq):
    (child_key1, child_branch1, branch_length1),(child_key1, child_branch1, branch_length1)=find_rooted_nodes(tree)

if __name__=='__main__':
    
    #print reduce_covariance(identity(10), 5)
    
    #print file_to_emp_cov('out_stem.cov',4)
    
    #from sys import exit
    
    #exit()
    from generate_prior_trees import generate_phylogeny
    from Rcatalogue_of_trees import *
    from Rtree_operations import create_trivial_tree, scale_tree
    tree2=scale_tree(generate_phylogeny(5,1),0.01)
    print pretty_string(tree2)
    print supplementary_text_ms_string()
    print tree_to_ms_command(tree2, 50)
    #print call_ms_string(supplementary_text_ms_string(), 'supp.txt')
    print call_ms_string(tree_to_ms_command(tree2, 50,50), 'tmp.txt')
    #cov= ms_to_treemix2('supp.txt', 20, 20,400)
    cov= ms_to_treemix2('tmp.txt', 50, 5,50)
    #cov2=calculate_covariance_matrix('tmp.txt', 50, 5,20)
    print cov
    #print cov2
    print make_covariance(tree2)
    print reduce_covariance(cov, 0)
    #print reduce_covariance(cov2, 0)
    print reduce_covariance(make_covariance(tree2),0)
    