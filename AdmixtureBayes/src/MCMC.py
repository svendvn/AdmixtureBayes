from Proposal_Function import prop
from scipy.stats import poisson
from likelihood import likelihood
from numpy.random import random
from math import exp
from tree_operations import get_number_of_admixes


def initialize_posterior(emp_cov):
    def posterior(x,pks={}):
        no_admixes=get_number_of_admixes(x)
#         tot_branch_length=0
        for branch in x:
            times=branch[2::2]
            less_than_zero=[1 for t in times if t<0]
            if sum(less_than_zero)>0:
                return -float('inf')
#         #print tot_branch_length
        prior=poisson.logpmf(no_admixs,0.1)
        likelihood=likelihood(x, emp_cov)
        pks['prior']=prior
        pks['likelihood']=likelihood
        return prior+likelihood
    return posterior
     

def one_jump(x, post, temperature, posterior_function,pks={}):
    
    newx,g1,g2,Jh,j1,j2=prop(x,pks)
    pks['prop_x']=newx
    pks['g1']=g1
    pks['g2']=g2
    pks['Jh']=Jh
    pks['j1']=j1
    pks['j2']=j2
    
#     print "x",x
#     print "newx",newx
#     print "g1", g1
#     print "g2",g2
#     print "Jh", Jh
#     print "j1",j1
#     print "j2", j2
    
    #print "newx",newx
    post_new=posterior_function(newx,pks)
    pks['prop_pos']=post_new
    
#     print "post_ratio", (post_new/post)
    
    mhr=exp(post_new-post)**temperature*g2*j2/j1/g1*Jh
    pks['mhr']=mhr
    #print "mhr", mhr
    #print newx
    #print post,post_new
    
#     print "mhr",mhr
    
    u=random()
    pks['U']=u
    if u<mhr:
        return newx,post_new
    return x,post

def basic_chain(start_tree, nodes, summaries, posterior_function, N=10000, verbose=True, sample_verbose_scheme=None):
    tree=start_tree
    post=posterior_function(tree)
    
    iteration_summary=[]
    
    verbose_list=_initialize_verbose_list(verbose_list, verbose, summaries)
    
    proposal_knowledge_scraper={}
        
    for i in range(N):
        new_tree,new_post=one_jump(tree, post, 1.0, posterior_function, proposal_knowledge_scraper)
        iteration_summary.append(_calc_and_print_summaries(verbose_list,
                                                           summaries,
                                                           new_tree=new_tree,
                                                           new_post=new_post,
                                                           old_post=post,
                                                           old_tree=old_tree,
                                                           iteration_number=i,**proposal_knowledge_scraper))
        tree=new_tree
        post=new_post
        
    
        
        
def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    res=[]
    iteration=kwargs['iteration_number']
    for s in summaries:
        save_num,print_num=sample_verbose_scheme.get(s.name, default=(-1,-1))
        save_bool = (save_num!=0) and (iteration % save_num==0) 
        print_bool = (print_num!=0) and (iteration % print_num==0)
        if save_bool or print_bool:
            val=s
    
def _initialize_verbose_list(verbose_list,verbose,summaries):
    if verbose_list is None and verbose:
        verbose_list=[]
        for n in range(len(summaries)):
            verbose_list.append((n,1))
    if not verbose:
        return []
    return verbose_list
    
    

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
    run_this(4)
    #cProfile.run('run_this(4)')
    
        

    