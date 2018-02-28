from Rtree_to_covariance_matrix import make_covariance
from Rtree_operations import (find_rooted_nodes, get_number_of_leaves, get_real_parents, pretty_string, get_no_leaves, 
                              node_is_non_admixture, node_is_admixture, node_is_leaf_node, node_is_coalescence, 
                              get_real_children_root, get_trivial_nodes, scale_tree_copy, get_leaf_keys,
                              time_adjust_tree, get_max_timing, get_all_branch_lengths, get_number_of_admixes)
from tree_statistics import get_timing, identifier_to_tree_clean, unique_identifier_and_branch_lengths
from load_data import read_data

from numpy import loadtxt, cov, array, mean, vstack, sum, identity, insert, hstack, vsplit, amin, sqrt, zeros, delete, ix_, ones
from numpy.linalg import det
from copy import deepcopy
from operator import itemgetter
import subprocess
from scipy.stats import wishart

def read_freqs(new_filename, locus_filter):
    with open(new_filename, 'r') as f:
        names=f.readline().split()
        allele_counts=[]
        pop_sizes=[]
        minors=[]
        total_sum=0
        for n,r in enumerate(f.readlines()):
            minor_majors=r.split()
            freqs=[]
            pop_sizes_SNP=[]
            for minor_major in minor_majors:
                minor, major= map(float,minor_major.split(','))
                total_sum+=major+minor
                freqs.append(float(minor)/float(major+minor))
                minors.append(minor)
                pop_sizes_SNP.append(major+minor)
            if locus_filter(freqs,pop_sizes, names):
                pop_sizes.append(pop_sizes_SNP)
                allele_counts.append(freqs)
    return allele_counts, names, pop_sizes, minors, total_sum

def make_uncompressed_copy(filename):
    take_copy_args=['cp', filename, filename+".tmp"]
    move_back_args=['mv', filename+'.tmp', filename]
    args=['gunzip', '-f', filename]
    new_filename='.'.join(filename.split('.')[:-1])
    subprocess.call(take_copy_args)
    subprocess.call(args)
    subprocess.call(move_back_args)
    return new_filename



def get_xs_and_ns_from_freqs(ps, npop, locus_filter):
    #npop=kwargs['sample_per_pop']
    mat=[]
    names=[]
    for k,p in ps.items():
        mat.append(p)
        names.append(k)
    mat=array(mat)
    ns=npop*ones(mat.shape)
    mat, ns= locus_filter.apply_filter(mat, ns, names)
    xs=mat*npop
    return xs,ns,names
    
def get_xs_and_ns_from_treemix_file(snp_file, locus_filter):
    if snp_file.endswith('.gz'):
        new_filename=make_uncompressed_copy(snp_file)
    else:
        new_filename=snp_file
    allele_freqs, names, ns, minors, total_sum= read_freqs(new_filename, locus_filter)
    allele_freqs=array(allele_freqs).T
    ns=array(ns).T
    xs=ns*allele_freqs
    return xs,ns,names

def order_covariance(xnn_tuple, outgroup=''):
    if not outgroup:
        return xnn_tuple
    xs,ns,names=xnn_tuple
    #print self.full_nodes, self.outgroup
    n_outgroup=next((n for n, e in enumerate(names) if e==outgroup))
    xs_o=xs[n_outgroup,:]
    ns_o=ns[n_outgroup,:]
    names_o=names[n_outgroup]
    xs=delete(xs, n_outgroup,0)
    ns=delete(ns, n_outgroup,0)
    names.remove(names_o)
    xs=insert(xs, 0, xs_o, axis=0)
    ns=insert(ns, 0, ns_o, axis=0)
    names=[names_o]+names
    return xs,ns,names
    
def _get_permutation(actual, target):
    val_to_index={val:key for key, val in enumerate(actual)}
    indices_to_get_target=[val_to_index[val] for val in target]
    return indices_to_get_target

def reorder_covariance(cov, names, full_nodes):
    indices=_get_permutation(names, full_nodes)
    return cov[ix_(indices,indices)]
    
def reorder_reduced_covariance(cov, names, full_nodes, outgroup=''):
    n1=len(names)
    n2=len(full_nodes)
    m=cov.shape[0]
    print 'cov', cov
    print 'names',names
    print 'full_nodes',full_nodes
    print outgroup
    assert n1==n2==(m+1), 'Unexpected input'
    names2=deepcopy(names)
    full_nodes2=deepcopy(full_nodes)
    names2.remove(outgroup)
    full_nodes2.remove(outgroup)
    indices=_get_permutation(names2, full_nodes2)
    print cov[ix_(indices, indices)]
    return cov[ix_(indices, indices)]

def tree_to_data_perfect_model(tree, df):
    m=make_covariance(tree)
    r=m.shape[0]
    m=wishart.rvs(df=r*df-1, scale=m/(r*df))
    return m

def normalise(m):
    return m-max(0,amin(m))

def file_to_emp_cov(filename, reduce_column=None, nodes=None):
    dat=[]
    with open(filename, 'r') as f:
        actual_nodes=f.readline().rstrip().split(" ")
        for l in f.readlines():
            print l
            n=map(float, l.split()[1:])
            if len(n)>1:
                dat.append(n)
    m=array(dat)
    print m
    mapping={val:key for key, val in enumerate(actual_nodes)}
    print 'mapping', mapping
    print 'nodes', nodes
    if nodes is not None:
        new_order=[mapping[node] for node in nodes]
        print 'new_order', new_order
        print 'm.shape', m.shape
        m=m[:, new_order][new_order]
    if reduce_column is not None:
        m=reduce_covariance(m, reduce_column)
        #m=normalise(m)
    return m

def emp_cov_to_file(m, filename='emp_cov.txt', nodes=None):
    if nodes is None:
        n=m.shape[0]
        nodes=get_trivial_nodes(n)
    with open(filename, 'w') as f:
        f.write(' '.join(nodes)+'\n')
        for i, node in enumerate(nodes):
            f.write(node+ ' '+ ' '.join(map(str, m[i]))+'\n')
    print 'wrote matrix to file', filename

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

def time_adjusted_tree_to_ms_command(time_adjusted_tree, sample_per_pop=50, nreps=2, 
                       theta=0.4, sites=500000, recomb_rate=1,
                       leaf_keys=None, final_pop_size=100.0):
    
    tree=deepcopy(time_adjusted_tree)
    if recomb_rate is None:
        rec_part=' -s '+str(sites)
    else:
        rec_part=' -r '+str(recomb_rate)+ ' '+str(sites)
    n=get_no_leaves(tree)
    callstring='ms '+str(sample_per_pop*n)+' '+str(nreps)+' -t '+ str(theta)+' ' +rec_part + ' '
    callstring+=' -I '+str(n)+' '+' '.join([str(sample_per_pop) for _ in xrange(n)])+' '
    times=get_max_timing(tree)
    print times
    tree=extend_branch_lengths(tree,times)
    print pretty_string(tree)
    if leaf_keys is None:
        leaf_keys= get_leaf_keys(tree)
    callstring+=construct_ej_es_string(tree, times, leaf_keys=leaf_keys, final_pop_size=final_pop_size)
    return callstring

def tree_to_ms_command(rtree, sample_per_pop=50, nreps=2, 
                       theta=1.0, sites=500000, recomb_rate=1,
                       leaf_keys=None, final_pop_size=100.0):
    tree=deepcopy(rtree)
    drift_sum=sum(get_all_branch_lengths(tree))
    if recomb_rate is None:
        rec_part=' -s '+str(sites)
    else:
        rec_part=' -r '+str(recomb_rate)+ ' '+str(sites)
    n=get_no_leaves(tree)
    callstring='ms '+str(sample_per_pop*n)+' '+str(nreps)+' -t '+ str(theta)+' ' +rec_part + ' '
    callstring+=' -I '+str(n)+' '+' '.join([str(sample_per_pop) for _ in xrange(n)])+' '
    times=get_timing(tree)
    #print times
    tree=extend_branch_lengths(tree,times)
    tuple_branch_lengths=get_all_branch_lengths(tree)
    count_sum=sum((x[1] for x in tuple_branch_lengths))
    tree=scaled_tupled_branches(tree, drift_sum/count_sum)
    times={k:v*drift_sum/count_sum for k,v in times.items()}
    #print pretty_string(tree)
    if leaf_keys is None:
        leaf_keys= get_leaf_keys(tree)
    callstring+=construct_ej_en_es_string(tree, times, leaf_keys=leaf_keys, final_pop_size=final_pop_size)
    
    
    #print tree
    #popsizes=[[calculate_pop_size(node[3])] if node_is_non_admixture(node) else [calculate_pop_size(node[3]), calculate_pop_size(node[4])] for key,node in tree.items()]    
    #pops=[p for l in popsizes for p in l]
    
    return callstring#,(min(pops),max(pops), max(times.values()))  #TO CHANGE BACK

def scaled_tupled_branches(tree, d):
    '''
    g
    '''
    for key in tree.keys():
        tree[key][3]=(tree[key][3][0], tree[key][3][1]*d)
        if tree[key][4] is not None:
            tree[key][4]=(tree[key][4][0], tree[key][4][1]*d)
    return tree
    
    
def extend_branch_lengths(tree, times):
    for key, node in tree.items():
        for n,parent_key in enumerate(get_real_parents(node)):
            pseudo_time=times[parent_key]-times[key]
            tree[key][3+n]=(tree[key][3+n], pseudo_time)
    return tree

def construct_ej_es_string(tree, times, leaf_keys, final_pop_size=1.0):
    s_times=sorted([(v,k) for k,v in times.items()])
    dic_of_lineages={(key,0):(n+1) for n,key in enumerate(leaf_keys)}
    print dic_of_lineages
    population_count=len(dic_of_lineages)
    res_string=''
    for time,key in s_times:
        if key=='r':
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            res_string+='-en '+str(time)+' '+str(dic_of_lineages[(key,0)])+' '+str(final_pop_size)+' '
            break
        node=tree[key]
        if node_is_coalescence(node):
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
        if node_is_admixture(node):
            population_count+=1
            i=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))[0]
            res_string+='-es '+str(time)+' '+str(i)+' '+str(1.0-node[2])+' '
            dic_of_lineages[(key,0)]=i
            dic_of_lineages[(key,1)]=population_count
    return res_string


def construct_ej_en_es_string(tree, times, leaf_keys, final_pop_size=1.0):
    s_times=sorted([(v,k) for k,v in times.items()])
    dic_of_lineages={(key,0):(n+1) for n,key in enumerate(leaf_keys)}
    #print dic_of_lineages
    population_count=len(dic_of_lineages)
    res_string=''
    for time,key in s_times:
        if key=='r':
            i,j=get_affected_populations(dic_of_lineages, get_real_children_root(tree, key))
            res_string+='-ej '+str(time)+' '+str(i)+' '+str(j)+' '
            dic_of_lineages[(key,0)]=j
            pop_size=final_pop_size
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
    return actual/drift*2

def call_ms_string(ms_string, sequence_file):
    with open(sequence_file, 'w') as f:
        print 'running', ms_string
        p = subprocess.Popen(ms_string, stdout=subprocess.PIPE, shell=True)
        line_number = 0
        for line in p.stdout.readlines():
            line_number += 1
            if line_number >= 5 and line and (line[0]=='0' or line[0]=='1'):
                #print len(line)
                #print line_number, line[:4]
                f.write(line)
            #else:
                #print line_number,':', line.rstrip()
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

    
    
    
    

def calculate_covariance_matrix2(file='tmp.txt', 
                                 samples_per_pop=20, 
                                 no_pops=4, 
                                 n_reps=1, 
                                 outgroup_number=None, 
                                 variance_correction=False,
                                 outgroup=True):
    data=[]
    with open(file, 'r') as f:
        for r in f.readlines():
            #print r[:5]
            data.append(map(int,list(r.rstrip())))
    new_data=[data[i:(i+samples_per_pop*no_pops)] for i in range(0,len(data),samples_per_pop*no_pops)]
    m=map(array, new_data)
    m=hstack(new_data)
    #print m.shape
    #alpha=0.05
    if outgroup:
        subtract_alleles=mean(m[outgroup_number*samples_per_pop:(outgroup_number+1)*samples_per_pop,:],axis=0)
    else:
        subtract_alleles=mean(m,axis=0)
    #other_alleles=mean(m[:outgroup_number*samples_per_pop], axis=0)/2+mean(m[(outgroup_number+1)*samples_per_pop:], axis=0)/2
    #sroot_homozygocity=sqrt(outgroup_alleles*(1.0-outgroup_alleles))
    #indices_of_positivity=[i for i,(s,o) in enumerate(zip(outgroup_alleles,other_alleles)) if s>alpha and s<(1.0-alpha) and o>alpha and o<(1.0-alpha)]
    #thinned_outgroup_alleles=outgroup_alleles[indices_of_positivity]
    #thinned_sroot_homozygocity=sroot_homozygocity[indices_of_positivity]
    #print 'SNPs', len(indices_of_positivity)
    ps=tuple([(mean(m[(i*samples_per_pop):((i+1)*samples_per_pop),  ], axis=0)-subtract_alleles) for i in xrange(no_pops)])
    p=vstack(ps)
    e_cov=cov(p)
    e_cov2=p.dot(p.T)/p.shape[1]
    #print 'e_cov2', e_cov2
    #print 'e_cov', e_cov
    e_cov=reduce_covariance(e_cov2, outgroup_number)
    return e_cov

def reduce_covariance(covmat, subtracted_population_index):
    reducer=insert(identity(covmat.shape[0]-1), subtracted_population_index, -1, axis=1)
    return reducer.dot(covmat).dot(reducer.T)

def ms_to_treemix2(filename='tmp.txt', samples_per_pop=20, no_pops=4, n_reps=1, filename2='tmp.treemix_in', treemix_files='tmp'):
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
    subprocess.call(['gzip','-f', filename2])
    return read_data(filename2_gz, blocksize=1000 ,outgroup='s1', noss=False, nodes=get_trivial_nodes(no_pops), outfile=treemix_files)

def treemix_to_cov(filename='treemix_in.txt.gz', outfile='not_used', reduce_also=False, reducer='', noss='not implemented', blocksize='unused', nodes=None):
    
    #unzip
    take_copy_args=['cp', filename, filename+".tmp"]
    move_back_args=['mv', filename+'.tmp', filename]
    args=['gunzip', '-f', filename]
    new_filename='.'.join(filename.split('.')[:-1])
    subprocess.call(take_copy_args)
    subprocess.call(args)
    subprocess.call(move_back_args)
    
    with open(new_filename, 'r') as f:
        names=f.readline().split()
        allele_counts=[]
        for r in f.readlines():
            minor_majors=r.split()
            freqs=[]
            for minor_major in minor_majors:
                minor, major= map(float,minor_major.split(','))
                freqs.append(float(minor)/float(major+minor))
            allele_counts.append(freqs)
    p=array(allele_counts)
    p=p.T
    
    mapping={val:key for key, val in enumerate(names)}
    if nodes is None:
        nodes=names
    new_order=[mapping[node] for node in nodes]
    p=p[new_order,:]
    
    if reduce_also:
        n_outgroup=next((n for n, e in enumerate(names) if e==reducer))
        print n_outgroup
        p=p-p[n_outgroup,:]
    m=p.dot(p.T)/p.shape[1]
    #print m
    
    
    if reduce_also:
        m=reduce_covariance(m, n_outgroup)
    return m
    
        
    


def ms_to_treemix3(filename='tmp.txt', samples_per_pop=20, no_pops=4, n_reps=1, filename2='tmp.treemix_in', nodes=None, 
                   convert_to_gz=True):
    if nodes is None:
        nodes=get_trivial_nodes(no_pops)
    total_sum=0
    total_number_of_genes=0
    with open(filename, 'r') as f:
        with open(filename2, 'w') as e:
            e.write(' '.join(nodes)+'\n')
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
                    total_sum+=sum(s_vec)
                    total_number_of_genes+=len(s_vec)*samples_per_pop
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
    muhat=float(total_sum)/float(total_number_of_genes)
    print 'muhat', muhat
    if not convert_to_gz:
        return filename2
    filename2_gz=filename2+'.gz'
    print filename2
    args=[['cp', '-T', filename2, filename2+'.tmp'],['gzip','-f', filename2], ['mv', filename2+'.tmp', filename2]]
    for arg in args:
        subprocess.call(arg)
    return filename2_gz
    
def ms_to_treemix(filename='tmp.txt', samples_per_pop=20, no_pops=4, n_reps=1, filename2='tmp.treemix_in', treemix_files='tmp'):
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
    subprocess.call(['gzip','-f', filename2])
    return read_data(filename2_gz, blocksize=10000 ,outgroup='s3', noss=True, outfile=treemix_files)

def trace_from_root(tree, init_freq):
    (child_key1, child_branch1, branch_length1),(child_key1, child_branch1, branch_length1)=find_rooted_nodes(tree)

def empirical_covariance_from_tree(tree, scaling=0.01, pop_size=20, reps=400):
    pass

def get_empirical_matrix(stree, factor=1.0, pop_size=20, reps=400):   
    tree= identifier_to_tree_clean(stree)
    ms_command=tree_to_ms_command(scale_tree_copy(tree, factor), pop_size, reps)
    print ms_command
    call_ms_string(ms_command, 'tmp.txt')
    empirical_covariance=ms_to_treemix2(filename='tmp.txt', samples_per_pop=pop_size, no_pops=get_number_of_leaves(tree), n_reps=reps, filename2='tmp.treemix_in')
    return reduce_covariance(empirical_covariance,0)

if __name__=='__main__':
    
    #print reduce_covariance(identity(10), 5)
    
    #print file_to_emp_cov('out_stem.cov',4)
    
    #from sys import exit
    
    #exit()
    from generate_prior_trees import generate_phylogeny
    from Rcatalogue_of_trees import *
    from Rtree_operations import create_trivial_tree, scale_tree
    tree2=scale_tree(generate_phylogeny(5,1),0.05)
    print pretty_string(tree2)
    print pretty_string(identifier_to_tree_clean(unique_identifier_and_branch_lengths(tree2)))
    print supplementary_text_ms_string()
    tree_good=generate_phylogeny(7)
    a=tree_to_ms_command(tree_good, 50,20)
    #print call_ms_string(a, 'supp.txt')
    b=time_adjusted_tree_to_ms_command(tree_good,50,20)
    #print call_ms_string(b, 'supp2.txt')
    #print call_ms_string(tree_to_ms_command(tree2, 50,20), 'tmp.txt')
    #cov= ms_to_treemix2('supp.txt', 20, 20,400)
    #cov= ms_to_treemix2('tmp.txt', 50, 5,20)
    #cov2=calculate_covariance_matrix('tmp.txt', 50, 5,20)
    #print cov
    #print cov2
    #print make_covariance(tree2)
    #print reduce_covariance(cov, 0)
    #print reduce_covariance(cov2, 0)
    #print reduce_covariance(make_covariance(tree2),0)
    
#     call_ms_string(a, 'sletmig/_ms_file.txt')
#     b1=calculate_covariance_matrix2('sletmig/_ms_file.txt', 50,7,20,0)
#     ms_to_treemix3('sletmig/_ms_file.txt', 50,7,20, filename2='sletmig/tmp_treemix_in.txt')
#     b2=read_data('sletmig/tmp_treemix_in.txt.gz', blocksize=1, nodes=None, noss=False, normalize=False, reduce_also=True, reducer='s1', return_muhat=False, outfile='tmp')
#     b3=treemix_to_cov('sletmig/tmp_treemix_in.txt.gz', 'outfile', True, 's1', noss=False, blocksize=None)
# 
# 
#     
#     c3=calculate_covariance_matrix2('sletmig/_ms.txt', 20,9,100,8)
#     
#     
#     from numpy import set_printoptions
#     set_printoptions(precision=4, suppress=True)
#     from load_data import read_data
#     cs=[read_data('sletmig/_treemix_in86.txt.gz', blocksize=1, nodes=None, noss=False, normalize=False, reduce_also=True, reducer=i, return_muhat=False, outfile='tmp') for i in range(9)]
#     c2=treemix_to_cov('sletmig/_treemix_in86.txt.gz', 'outfile', True, 'out', noss=False, blocksize=None)
#     
#     print b2
#     print b3
#     print b1
#     
#     print c2
#     print c3
#     for c in cs:
#         print c
        
    testdata=['0010','\n','0011','\n','0000','\n','1110','\n','1000','\n','1100','\n','1100','\n','1111','\n','0111']
    snp_file=['0,3 3,0 2,1','\n','0,3 2,1 3,0', '\n', '2,1 1,2 2,1','\n','1,2 0,3 2,1']
    with open('sletmig/ms.txt', 'w') as f:
        for a in testdata:
            f.write(a)
    with open('sletmig/treemix_in.txt' , 'w') as f:
        f.write('s1 s2 out'+'\n')
        for a in snp_file:
            f.write(a)
    subprocess.call(['gzip', '-f','-k', 'sletmig/treemix_in.txt'])
    b1=calculate_covariance_matrix2('sletmig/ms.txt', 3,3,1,2)
    b2=[read_data('sletmig/treemix_in.txt.gz', blocksize=4, nodes=['out','s1','s2'], noss=True, normalize=False, reduce_also=True, reducer=i, return_muhat=False, outfile='tmp') for i in range(3)]
    b3=treemix_to_cov('sletmig/treemix_in.txt.gz', 'outfile', True, 'out', noss=False, blocksize=None)

    for b in b2:
        print b
    print b3
    print b1
    
    