from tree_statistics import identifier_to_tree_clean
from tree_plotting import plot_as_directed_graph, pretty_string
#from sphinx.util.nodes import _new_copy


true_tree_s= 'w.w.w.w.w.w.a.a.w-c.w.c.c.w.c.5.0.w.3.2-c.w.w.0.c.4.w-c.w.0.c.4-w.c.1-c.0;0.07-0.974-1.016-0.089-0.81-0.086-1.499-0.052-1.199-2.86-0.403-0.468-0.469-1.348-1.302-1.832-0.288-0.18-0.45-0.922-2.925-3.403;0.388-0.485'
true_tree=identifier_to_tree_clean(true_tree_s)

wrong_trees_s=['w.w.a.w.w.a.a.a.w-c.w.c.c.w.w.c.0.w.w.6.3.2-c.w.w.0.w.c.5.w.w-c.w.0.c.4.w.w-c.w.c.4.0-w.c.1-c.0;0.828-0.21-0.197-0.247-0.568-1.06-0.799-1.162-2.632-2.001-0.45-1.048-0.834-0.469-0.191-2.759-0.871-1.896-0.473-0.019-1.236-0.287-0.179-0.981-0.456-0.91-2.114-3.368;0.655-0.506-0.389-0.23',
               'w.w.w.w.w.w.a.a.w-w.w.w.c.w.c.5.a.w.3.a-c.w.c.w.c.4.w.w.0.2.a-w.w.w.w.c.c.4.5.w-c.c.w.w.1.0.w-c.w.w.0.w-c.w.0.w-a.w.w-c.w.0.w-c.w.0-c.0;0.387-0.087-0.806-0.082-2.062-0.803-0.122-0.544-0.061-0.733-0.474-1.342-0.871-0.798-0.753-0.288-0.024-0.174-0.754-0.282-0.45-0.924-0.416-1.081-0.467-1.296-1.171-0.54-1.944-0.258-8.813-0.76-0.073-3.416;0.388-0.467-0.098-0.185-0.019-0.44']
wrong_trees=[identifier_to_tree_clean(tree) for tree in wrong_trees_s]

plot_as_directed_graph(true_tree,  drawing_name= 'tmp0.bmp')
plot_as_directed_graph(wrong_trees[0], drawing_name = 'tmp1.bmp')
print pretty_string(wrong_trees[0])
t=wrong_trees[0]

from Rproposal_admix import deladmix

pks={}
from Rtree_to_covariance_matrix import make_covariance
from posterior import initialize_big_posterior

true_cov=make_covariance(true_tree)
posterior_f=initialize_big_posterior(true_cov, M=10000)
nt, f,b=deladmix(t,pks=pks, fixed_remove=('a1',1))
plot_as_directed_graph(nt)

new_likelihood_value, new_prior_value, (new_branch_prior, new_no_admix_prior, new_admix_prop_prior, new_top_prior), new_covariance= posterior_f((nt,0))
old_likelihood_value, old_prior_value, (old_branch_prior, old_no_admix_prior, old_admix_prop_prior, old_top_prior), old_covariance= posterior_f((t,0))

print new_likelihood_value, old_likelihood_value, new_likelihood_value-old_likelihood_value
from numpy import get_printoptions, set_printoptions
set_printoptions(precision=3, suppress=True)
print new_covariance
print old_covariance
print new_covariance-old_covariance
print true_cov
print true_cov-old_covariance
print true_cov-new_covariance
