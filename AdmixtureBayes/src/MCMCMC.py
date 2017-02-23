from MCMC import basic_chain, initialize_posterior
from pathos.multiprocessing import Pool

def basic_chain_unpacker(args):
    """Convert `f([1,2])` to `f(1,2)` call.
    Thank you stackexchange"""
    return basic_chain(*args)

def MCMCMC():
    pool = Pool()
    a_args = [1,2,3]
    second_arg = 1
    pool.map(basic_chain, itertools.izip(a_args, itertools.repeat(second_arg)))
    
from scipy.stats import poisson
from likelihood import likelihood
from numpy.random import random
from math import exp
from tree_operations import get_number_of_admixes
from summary import *
import itertools
    
    
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

pool = Pool()
a_args = [x,summaries, posterior_function, 1000, sample_verbose_scheme, 10, 0]
temps=range(1,4)
pool.map(basic_chain_unpacker, itertools.izip(itertools.repeat(a_args), temps))