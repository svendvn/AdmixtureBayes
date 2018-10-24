AdmixtureBayes run --input_file allele_counts.txt --outgroup s4 --bootstrap_blocksize 2
AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt 
AdmixtureBayes plot --posterior_distribution_file posterior_distributions.csv --plot consensus_trees --write_ranking_to_file rankings.txt
