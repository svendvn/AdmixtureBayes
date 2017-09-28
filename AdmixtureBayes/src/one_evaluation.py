
def one_evaluation(starting_tree, 
                   posterior_function,
                   result_file):
    likelihood_val, prior_val= posterior_function(starting_tree)
    posterior_val = likelihood_val+prior_val
    
    with open(result_file, 'w') as f:
        f.write(' '.join([likelihood_val, prior_val, posterior_val]))
    