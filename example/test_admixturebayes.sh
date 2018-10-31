AdmixtureBayes run --input_file big_allele_counts.txt --outgroup out --prefix svend
AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt 
AdmixtureBayes plot --posterior_distribution_file posterior_distributions.csv --plot consensus_trees --write_ranking_to_file rankings.txt
