AdmixtureBayes run --input_file big_allele_counts.txt --outgroup out
AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt --outgroup_name out 
AdmixtureBayes plot --posterior_distribution_file posteriors_small.csv --plot consensus_trees --write_ranking_to_file rankings.txt
#AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt --outgroup_name out --subnodes s1 s2 s3 --result_file posteriors_small.csv
#AdmixtureBayes qp --posterior_distribution_file posterior_distributions.csv 
