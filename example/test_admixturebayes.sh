AdmixtureBayes run --input_file big_allele_counts.txt --outgroup out --n 500
AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt  --total 300 
AdmixtureBayes plot --posterior_distribution_file posterior_distributions.csv --plot consensus_trees --write_ranking_to_file rankings.txt
AdmixtureBayes plot --posterior_distribution_file posterior_distributions.csv --plot top_node_trees --write_ranking_to_file nrankings.txt
AdmixtureBayes plot --posterior_distribution_file posterior_distributions.csv --plot top_trees --write_ranking_to_file trankings.txt
AdmixtureBayes plot --posterior_distribution_file posterior_distributions.csv --plot estimates

AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt --outgroup_name out --subnodes s1 s2 s3 --result_file posterior_distributions.csv
AdmixtureBayes qp --posterior_distribution_file posterior_distributions.csv 

mkdir -p variations

#AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt \
#--result_file variations/posteriors_small_no_outgroup.csv --total 300 
#AdmixtureBayes plot --posterior_distribution_file variations/posteriors_small_no_outgroup.csv --plot top_trees \
#--write_ranking_to_file trankings.txt --popup --prefix variations/no_outgroup

#AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt \
#--result_file variations/posteriors_small_no_outgroup_subnodes.csv --total 300 --subnodes s1 s2 s3 
#AdmixtureBayes plot --posterior_distribution_file variations/posteriors_small_no_outgroup_subnodes.csv --plot top_trees \
#--write_ranking_to_file trankings.txt --popup --prefix variations/no_outgroup_subnodes


#AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt \
#--result_file variations/posteriors_small.csv --total 300  --outgroup_name out --subnodes s1 s2 s3
#AdmixtureBayes plot --posterior_distribution_file variations/posteriors_small.csv --plot top_trees \
#--write_ranking_to_file trankings.txt --popup --prefix variations/subnodes_wo_outgroup


#AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt \
#--result_file variations/posteriors_small.csv --total 300   --outgroup_name out --subnodes s1 s2 out
#AdmixtureBayes plot --posterior_distribution_file variations/posteriors_small.csv --plot top_trees \
#--write_ranking_to_file trankings.txt --popup --prefix variations/subnodes_with_outgroup

