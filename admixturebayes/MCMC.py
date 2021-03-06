from Proposal_Function import prop_flat
from scipy.stats import poisson, norm, multivariate_normal
from likelihood import likelihood
from numpy.random import random
from math import exp, log, fabs
from tree_operations import get_number_of_admixes, illegal_branch_length
from summary import *
from multiprocessing import Queue, Process
from prior import prior
from tree_warner import check
from Rtree_operations import pretty_string, scale_tree_copy



            
def one_jump(x, post, temperature, posterior_function, proposal, pks={}):
    check_old = False
    if check_old:
        # pass
        post_old = posterior_function(x)

        if fabs(post_old[1] - post[1]) > 0.00001 :
            print 'The likelihood to the old tree is not correct: ' + str(fabs(post_old[1] - post[1]))
            print posterior_function.base_r
    
    newx,g1,g2,Jh,j1,j2=proposal(x,pks)
    pks['proposed_tree']=newx[0]
    pks['g1']=g1
    pks['g2']=g2
    pks['Jh']=Jh
    pks['j1']=j1
    pks['j2']=j2
    
    post_new=posterior_function(newx,pks)
    pks['proposed_posterior']=post_new


    #print temperature
    likelihood_old, prior_old = post[:2]
    likelihood_new, prior_new = post_new[:2]
    #print 'Posterior Matrix:'
    #print likelihood_old, prior_old
    #print likelihood_new, prior_new
    
    #for key,val in pks.items():
    #    print key, '=', val
    
    if g2<=0 or j2<=0:
        logmhr=-float('inf')
    else:
        logmhr=(likelihood_new-likelihood_old)/temperature+(prior_new-prior_old)+log(g2)+log(j2)-log(j1)-log(g1)+log(Jh)
    if logmhr>100:
        mhr=float('inf')
    else:
        mhr=exp(logmhr)
        
    #print post_new, post, post_new-post, exp(post_new-post), mhr
    pks['mhr']=mhr
    
    u=random()
    pks['U']=u
    proposal.adapt(mhr, u, post_new, post, temperature)
    if u<mhr:
        return newx,post_new
    return x,post


def basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=10000, 
                sample_verbose_scheme=None, overall_thinning=1, i_start_from=0, 
                temperature=1.0, proposal_update=None, multiplier=None, check_trees=False, 
                appending_result_file=None, appending_result_frequency=10):
    if proposal_update is not None:
        proposal.wear_exportable_state(proposal_update)
    
    x=start_x
    if check_trees:
        check(x[0])
    if post is None:
        post=posterior_function(x)
    
    iteration_summary=[]
    #print 'random', random()
    count=0
    from_count=0
    
    if appending_result_file is not None:
        with open(appending_result_file, 'w') as f:
            f.write(",".join(['iteration'] + [s.name for s in summaries])+'\n')
        
    for i in range(i_start_from,i_start_from+N):
        proposal_knowledge_scraper={}
        new_x,new_post=one_jump(x, post, temperature, posterior_function, proposal, proposal_knowledge_scraper)
        if overall_thinning!=0 and i%overall_thinning==0:
            if multiplier is None:
                iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               tree=new_x[0],
                                                               add=new_x[1],
                                                               posterior=new_post,
                                                               old_post=post,
                                                               old_tree=x[0],
                                                               iteration_number=i,**proposal_knowledge_scraper))
            else:
                iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               tree=scale_tree_copy(new_x[0], 1.0/multiplier),
                                                               add=new_x[1],
                                                               posterior=new_post,
                                                               old_post=post,
                                                               old_tree=scale_tree_copy(x[0],1.0/multiplier),
                                                               iteration_number=i,**proposal_knowledge_scraper))
            if appending_result_file is not None:
                count+=1
                if count % appending_result_frequency==0:
                    with open(appending_result_file, 'a') as f:
                        for n,params in enumerate(iteration_summary[from_count:]):
                            f.write(",".join(map(str, params))+'\n')
                    from_count=count
        x=new_x
        post=new_post
        if check_trees:
            #print pretty_string(tree)
            check(x[0], proposal_knowledge_scraper)
    
    #print iteration_summary
    #print zip(*iteration_summary)
    return x, post, zip(*iteration_summary),proposal.get_exportable_state()
        
        
        
def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    iteration=kwargs['iteration_number']
    res=[iteration]
    for s in summaries:
        save_num,print_num=sample_verbose_scheme.get(s.name, (0,0))
        save_bool = (save_num!=0) and (iteration % save_num==0) 
        print_bool = (print_num!=0) and (iteration % print_num==0)
        if save_bool or print_bool:
            val=s(**kwargs)
            if print_bool:
                print str(iteration)+'. '+ s.pretty_print(val)
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

    