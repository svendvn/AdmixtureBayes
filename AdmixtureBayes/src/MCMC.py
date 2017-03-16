from Proposal_Function import prop_flat
from scipy.stats import poisson
from likelihood import likelihood
from numpy.random import random
from math import exp
from tree_operations import get_number_of_admixes, illegal_branch_length
from summary import *
from multiprocessing import Queue, Process
from prior import prior


def initialize_posterior(emp_cov):
    def posterior(x,pks={}):
        #print tot_branch_length
        prior_val=prior(x)
        if prior_val==-float('inf'):
            return prior_val
        likelihood_value=likelihood(x, emp_cov)
        pks['prior']=prior_val
        pks['likelihood']=likelihood_value
        return prior_val+likelihood_value
    return posterior
     

def one_jump(x, post, temperature, posterior_function, proposal, pks={}):
    
    newx,g1,g2,Jh,j1,j2=proposal(x,pks)
    pks['proposed_tree']=newx
    pks['g1']=g1
    pks['g2']=g2
    pks['Jh']=Jh
    pks['j1']=j1
    pks['j2']=j2
    
    post_new=posterior_function(newx,pks)
    pks['proposed_posterior']=post_new
    
    mhr=exp(post_new-post)**temperature*g2*j2/j1/g1*Jh
    
    pks['mhr']=mhr
    
    u=random()
    pks['U']=u
    #proposal.update(mhr, u, post_new, post, temperature)
    if u<mhr:
        return newx,post_new
    return x,post


def basic_chain(start_tree, summaries, posterior_function, proposal, post=None, N=10000, sample_verbose_scheme=None, overall_thinning=1, i_start_from=0, temperature=1.0, proposal_update=None):
    if proposal_update is not None:
        proposal.update(proposal_update)
    
    tree=start_tree
    if post is None:
        post=posterior_function(tree)
    
    iteration_summary=[]
    
    proposal_knowledge_scraper={}
        
    for i in range(i_start_from,i_start_from+N):
        new_tree,new_post=one_jump(tree, post, temperature, posterior_function, proposal, proposal_knowledge_scraper)
        if overall_thinning!=0 and i%overall_thinning==0:
            iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               tree=new_tree,
                                                               posterior=new_post,
                                                               old_post=post,
                                                               old_tree=tree,
                                                               iteration_number=i,**proposal_knowledge_scraper))
        tree=new_tree
        post=new_post
    
    return tree, post, zip(*iteration_summary), None#,proposal.get_update()
        
        
        
def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    res=[]
    iteration=kwargs['iteration_number']
    for s in summaries:
        save_num,print_num=sample_verbose_scheme.get(s.name, (0,0))
        save_bool = (save_num!=0) and (iteration % save_num==0) 
        print_bool = (print_num!=0) and (iteration % print_num==0)
        if save_bool or print_bool:
            val=s(**kwargs)
            if print_bool:
                print s.pretty_print(val)
            if save_bool:
                res.append(val)
            else:
                res.append(None)
        else:
            res.append(None)
    return res
    
    

if __name__=="__main__":
    
    def run_this(N):
        
        from tree_operations import make_flat_list_no_admix
        
        from numpy import diag
        tree_flatter_list=make_flat_list_no_admix(N)
        nodes=["s"+str(i) for i in range(1,N+1)]
        emp_cov=diag([0.5]*N)
        emp_cov[2,1]=emp_cov[1,2]=0.2
    
        sigma=1
        x=tree_flatter_list
        posterior=initialize_posterior(emp_cov)
        post=posterior(x)
        
        #print "start post_true"
        #print posterior(true_tree)
        #print "end post_true"
        
    
        admixes=[]
        posteriors=[]
        for i in range(10000):
            x,post=one_jump(x,post, 1.0,posterior)
            posteriors.append(post)
            tot_length=sum(map(len, x))
            no_admixs=(tot_length-len(x)*3)/4
            admixes.append(no_admixs)
           # if i%1000==0:
                ##print i
        from collections import Counter
        print Counter(admixes)
        #print posteriors
    
    import cProfile
    #run_this(4)
    #cProfile.run('run_this(4)')
    from tree_operations import make_flat_list_no_admix
        
    from numpy import diag
    N=5
    tree_flatter_list=make_flat_list_no_admix(N)
    nodes=["s"+str(i) for i in range(1,N+1)]
    emp_cov=diag([0.5]*N)
    emp_cov[2,1]=emp_cov[1,2]=0.2
    
    x=tree_flatter_list
    posterior_function=initialize_posterior(emp_cov)
    summaries=[s_variable('posterior'), s_variable('mhr'), s_branch_length()]
    sample_verbose_scheme={'posterior':(1,1),
                           'branch_length':(10,0),
                           'mhr':(1,0)}
    #rd=basic_chain(x, summaries, posterior_function, N=10000, sample_verbose_scheme=sample_verbose_scheme, overall_thinning=10)
    cProfile.run('basic_chain(x, summaries, posterior_function, N=10000, sample_verbose_scheme=sample_verbose_scheme, overall_thinning=10)')

    