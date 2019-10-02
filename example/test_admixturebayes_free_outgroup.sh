#here we allow admixture events to the outgroup by writing outgroup_type Free. To increase chance of convergence we also set a higher MCMC_chains. Unrelated, we save the files in another folder
mkdir -p free_outgroup 
AdmixtureBayes run --input_file big_allele_counts.txt --outgroup out --outgroup_type Free --MCMC_chains 16 --prefix free_outgroup/_ --result_file free_outgroup/free_outgroup.csv
AdmixtureBayes posterior --input_file free_outgroup/free_outgroup.csv --covariance_matrix_file free_outgroup/_covariance_and_multiplier.txt --outgroup_name out --result_file free_outgroup/posteriors_small.csv --total 300 
AdmixtureBayes plot --posterior_distribution_file free_outgroup/posteriors_small.csv --plot consensus_trees --write_ranking_to_file rankings.txt --prefix free_outgroup/
AdmixtureBayes plot --posterior_distribution_file free_outgroup/posteriors_small.csv --plot top_node_trees --write_ranking_to_file nrankings.txt --prefix free_outgroup/
AdmixtureBayes plot --posterior_distribution_file free_outgroup/posteriors_small.csv --plot top_trees --write_ranking_to_file trankings.txt --prefix free_outgroup/
AdmixtureBayes plot --posterior_distribution_file free_outgroup/posteriors_small.csv --plot estimates --prefix free_outgroup/

#to make results more interpretable we set the rerooting using two different methods
mkdir -p free_outgroup/rerooted_force
mkdir -p free_outgroup/rerooted_ignore
AdmixtureBayes posterior --input_file free_outgroup/free_outgroup.csv --covariance_matrix_file free_outgroup/_covariance_and_multiplier.txt --outgroup_name out \
                         --result_file free_outgroup/rerooted_force/posteriors_small.csv --total 300 --reroot s4 --reroot_error force 
AdmixtureBayes plot --posterior_distribution_file free_outgroup/rerooted_force/posteriors_small.csv --plot top_trees --write_ranking_to_file trankings.txt --prefix free_outgroup/rerooted_force/
AdmixtureBayes posterior --input_file free_outgroup/free_outgroup.csv --covariance_matrix_file free_outgroup/_covariance_and_multiplier.txt --outgroup_name out \
                         --result_file free_outgroup/rerooted_ignore/posteriors_small.csv --total 300 --reroot s2 --reroot_error ignore
AdmixtureBayes plot --posterior_distribution_file free_outgroup/rerooted_ignore/posteriors_small.csv --plot top_trees --write_ranking_to_file trankings.txt --prefix free_outgroup/rerooted_ignore/

#combining it with subnodes
#mkdir -p free_outgroup/sub_with_outgroup
#mkdir -p free_outgroup/sub_without_outgroup
#AdmixtureBayes posterior --input_file free_outgroup/free_outgroup.csv --covariance_matrix_file free_outgroup/_covariance_and_multiplier.txt --outgroup_name out \
#                         --result_file free_outgroup/sub_with_outgroup/posteriors_small.csv --total 300 --reroot out --reroot_error force --subnodes s1 s2 out 
#AdmixtureBayes plot --posterior_distribution_file free_outgroup/sub_with_outgroup/posteriors_small.csv --plot top_trees --write_ranking_to_file trankings.txt --prefix free_outgroup/sub_with_outgroup/
#AdmixtureBayes posterior --input_file free_outgroup/free_outgroup.csv --covariance_matrix_file free_outgroup/_covariance_and_multiplier.txt --outgroup_name out \
#                         --result_file free_outgroup/sub_without_outgroup/posteriors_small.csv --total 300 --reroot out --reroot_error force --subnodes s1 s2 s3 
#AdmixtureBayes plot --posterior_distribution_file free_outgroup/sub_without_outgroup/posteriors_small.csv --plot top_trees --write_ranking_to_file trankings.txt --prefix free_outgroup/sub_without_outgroup/
#AdmixtureBayes posterior --input_file free_outgroup/free_outgroup.csv --covariance_matrix_file free_outgroup/_covariance_and_multiplier.txt --outgroup_name out \
#                         --result_file free_outgroup/sub_without_outgroup/posteriors_small.csv --total 300 --subnodes s1 s2 s3 --faster 
#AdmixtureBayes plot --posterior_distribution_file free_outgroup/sub_without_outgroup/posteriors_small.csv --plot top_trees --write_ranking_to_file trankings.txt --prefix free_outgroup/sub_without_outgroup/
#AdmixtureBayes plot --posterior_distribution_file free_outgroup/sub_without_outgroup/posteriors_small.csv --plot estimates --prefix free_outgroup/sub_without_outgroup/

#

#AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt --outgroup_name out --subnodes s1 s2 s3 --result_file posteriors_small.csv
#AdmixtureBayes qp --posterior_distribution_file posterior_distributions.csv 
