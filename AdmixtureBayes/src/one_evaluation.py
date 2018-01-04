
def one_evaluation(starting_tree, 
                   posterior_function,
                   result_file):
    likelihood_val, prior_val= posterior_function(starting_tree, verbose=True)
    posterior_val = likelihood_val+prior_val
    likelihood_val2= posterior_function.get_max_likelihood(verbose=True)
    likelihood_val3= posterior_function.get_non_empirical_max_likelihood(starting_tree, verbose=True)
    other_statistics= posterior_function.get_size_diff(starting_tree)
    
    with open(result_file, 'w') as f:
        f.write(' '.join(['starting_tree']+map(str,[likelihood_val, prior_val, posterior_val]))+'\n')
        f.write(' '.join(['max_likelihood']+map(str,[likelihood_val2, '.', '.']))+'\n')
        f.write(' '.join(['pmax_likelihood']+map(str,[likelihood_val3, '.', '.']))+'\n')
        for name, stat in zip(['median_ratio', 'relative_sign', 'frobenius'], other_statistics):
            f.write(' '.join([name]+map(str,[stat, '.', '.']))+'\n')
    