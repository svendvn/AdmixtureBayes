from tree_statistics import identifier_to_tree_clean,topological_identifier_to_tree_clean
import tree_generation_laboratory as tgl
from generate_prior_trees import generate_phylogeny
from tree_plotting import plot_as_directed_graph
import numpy as np
from Rtree_operations import get_leaf_keys
from Rtree_to_covariance_matrix import make_covariance
from scipy.stats import wishart as scip_wishart
from scipy.stats import chi2
import summary
from construct_proposal_choices import make_proposal
from posterior import posterior_class
from MCMC import basic_chain
import tree_statistics
import pandas as pd
from MCMCMC import MCMCMC
from temperature_scheme import fixed_geometrical
import collections
import tree_to_data 
from construct_filter_choices import make_filter
import Rtree_operations

R=np.array([[-1]*3, [1,0,0],[0,1,0],[0,0,1]]).T
nodes=['s1','s2','s3','s4']

def cov(p):
    return p.dot(p.T)/p.shape[1]

def scale_tree(tree, mult):
    return tree_statistics.unique_identifier_and_branch_lengths(Rtree_operations.scale_tree(identifier_to_tree_clean(tree),mult),nodes)

def run_ms_and_convert_to_treemix_format(ms_tree, samples_per_population=10, nreps=2, no_pops=4, run_ms=True, filename='ms.txt'):
    if run_ms:
        tree_to_data.call_ms_string(ms_tree, filename)
    filename2='snp_data.txt'
    filename2=tree_to_data.ms_to_treemix3(filename, 
                    samples_per_pop=samples_per_population, 
                    no_pops=no_pops, 
                    n_reps=nreps, 
                    filename2=filename2,
                    nodes=['s'+str(i) for i in range(1,no_pops+1)],
                    convert_to_gz=False)
    locus_filter=make_filter()
    xs,ns,names= tree_to_data.get_xs_and_ns_from_treemix_file(filename2, locus_filter)
    dat_frame={}
    for x,n,name in zip(np.array(xs), np.array(ns), names):
        dat_frame[name]=[str(int(xi))+','+str(int(ni-xi)) for xi,ni in zip(x,n)]
    return pd.DataFrame(dat_frame)

def convert_to_allele_frequency(snp_data):
    snp_data.to_csv('tmp.csv', sep=' ', index=False)
    locus_filter=make_filter()
    xs,ns,names=tree_to_data.get_xs_and_ns_from_treemix_file('tmp.csv', locus_filter)
    dat_frame={}
    for x,n,name in zip(np.array(xs), np.array(ns), names):
        dat_frame[name]=[float(xi)/ni for xi,ni in zip(x,n)]
    return pd.DataFrame(dat_frame)
    
def _cheat_and_load_previously_created_data(filename2):
    return run_ms_and_convert_to_treemix_format(None, run_ms=False, filename='ms_backup.txt')
    

def tree_to_ms_command(stree, samples_per_population=10, snps=1000000):
    nreps=snps//500000
    tree=identifier_to_tree_clean(stree)
    return tree_to_data.tree_to_ms_command(tree, sample_per_pop=samples_per_population, nreps=nreps)
    
def simulate_tree(*args):
    np.random.seed(3232)
    return tgl.simulate_tree(*args)

def visualize_topology(stree):
    numba=str(np.random.randint(0,1000000))
    filename='tree'+numba+'.png'
    if ';' in stree:
        plot_as_directed_graph(identifier_to_tree_clean(stree), drawing_name= filename, popup=False)
        
    else:
        plot_as_directed_graph(topological_identifier_to_tree_clean(stree), drawing_name= filename, popup=False)
    return filename

def tree_to_covariance(stree):
    tree=identifier_to_tree_clean(stree)
    nodes=sorted(get_leaf_keys(tree))
    return make_covariance(tree, node_keys=nodes)

def wishart(covariance,mean,df):
    if isinstance(covariance, int) or isinstance(covariance, float):
        return chi2.pdf(covariance, scale=mean/df, df=df)
    return scip_wishart.pdf(covariance, scale=mean/df, df=df)

def logwishart(covariance,mean,df):
    if isinstance(covariance, int) or isinstance(covariance, float):
        return chi2.logpdf(covariance, scale=mean/df, df=df)
    return scip_wishart.logpdf(covariance, scale=mean/df, df=df)

class options_object():
    
    def __init__(self,add=0, chains=1):
        self.deladmix=1
        self.addadmix=1
        self.rescale=1
        self.regraft=0
        self.rescale_add=int(add)
        self.rescale_admix=1
        self.rescale_admix_correction=0
        self.rescale_constrained=1
        self.rescale_marginally=0
        self.sliding_regraft=1
        self.sliding_rescale=0
        self.cancel_preserve_root_distance=False
        self.no_add=(not add)
        self.MCMC_chains=chains

def pretty_dataframe(result):
    return pd.DataFrame(map(list,zip(*result[:3])), columns=['iteration','posterior','tree'])        

def mcmcmc(observed_covariance, df , outgroup=False, chains=8, its=[50]*100):
    nodes=['s'+str(i+1) for i in range(observed_covariance.shape[0])]
    start_x=identifier_to_tree_clean(simulate_tree(4,0)),0
    summaries=[summary.s_posterior(), 
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_variable('add', output='double'),
               summary.s_no_admixes(),]
    options=options_object(outgroup, chains=chains)
    proposal=make_proposal(options)
    posterior_function=posterior_class(observed_covariance, M=df, nodes=nodes)
    sample_verbose_scheme=[{'posterior':(1,200), 'tree':(1,0),'add':(1,200),'no_admixes':(1,200)}]+[{s.name:(1,0) for s in summaries}]*(chains-1)

    res=MCMCMC(starting_trees=[(identifier_to_tree_clean(simulate_tree(4,0)),0) for _ in range(chains)], 
       posterior_function= posterior_function,
       summaries=summaries, 
       temperature_scheme=fixed_geometrical(800,chains), 
       printing_schemes=sample_verbose_scheme, 
       iteration_scheme=its, 
       overall_thinnings=40, 
       proposal_scheme= proposal, 
       cores=chains, 
       no_chains=chains,
       multiplier=None,
       result_file=None,
       store_permuts=False, 
       stop_criteria=None)
    res=res.loc[res.layer==0,['iteration','posterior','tree','no_admixes']]
    return res
      
def mcmc(observed_covariance, df, outgroup=False):
    nodes=['s'+str(i+1) for i in range(observed_covariance.shape[0])]
    start_x=identifier_to_tree_clean(simulate_tree(4,0)),0
    summaries=[summary.s_posterior(), 
               summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               summary.s_variable('add', output='double'),
               summary.s_no_admixes(),]
    options=options_object(outgroup)
    proposal=make_proposal(options)[0]
    posterior_function=posterior_class(observed_covariance, M=df, nodes=nodes)
    sample_verbose_scheme={'posterior':(1,200), 'tree':(1,0),'add':(1,200),'no_admixes':(1,200)}
    a=basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=5000, 
                sample_verbose_scheme=sample_verbose_scheme, overall_thinning=100, i_start_from=0, 
                temperature=1.0, proposal_update=None, multiplier=None, check_trees=False, 
                appending_result_file=None, appending_result_frequency=10)
    return a[2]

def top_topology_frequencies(result, top=3):
    trees=result['tree']
    topologies=[t.split(';')[0] for t in trees]
    ad_c=collections.Counter(topologies)
    tops=ad_c.most_common(top)
    return tops

def get_covariance_matrix_example_from_topology(result, topology):
    trees=result['tree']
    #print trees
    topologies=np.array([t.split(';')[0] for t in trees])
    #print topologies==topology
    relevant_trees=trees.loc[topologies==topology]
    #print relevant_trees
    return tree_to_covariance(relevant_trees.iloc[-1])
    
    
if __name__=='__main__':
    tree=simulate_tree(4,2)
    tree_covariance=tree_to_covariance(tree)
    r=mcmcmc(observed_covariance=tree_covariance,df=1000, its=[50]*800)    