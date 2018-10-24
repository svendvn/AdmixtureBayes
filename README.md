# AdmixtureBayes
## Purpose
AdmixtureBayes estimates posterior samples of admixture graphs (phylogenies incorporating admixture events) given an allele count file.

AdmixtureBayes is also shipped with commands to summarize and plot the results

## Installation


### Linux 64-bit

For Linux 64-bit, install the program in a python 2.7 environment using the command
```bash
$ pip install dist/AdmixtureBayes-0.1-cp27-cp27mu-linux_x86_64.whl
```
It will install the necessary python dependencies. Installation of the pygraphviz package is likely to fail with pip. Install it in another way, e.g. by conda
```bash
$ conda install -c anaconda pygraphviz 
```

### Other OS

For other than linux 64-bit, the only option (for now) is to install a version without the C-accelerated likelihood
```bash
$ pip install dist/AdmixtureBayes-0.1-cp27-cp27mu-linux_x86_64.whl
```
 
 
## Input file

The input for AdmixtureBayes is an allele count file in the exact same format as used by TreeMix.

```bash
$ cat allele_counts.txt
population1 population2 population3 population4
3,9 6,6 0,12 8,1
12,0 12,0 12,0 0,11
...
5,7 2,10 1,11 6,6
```

where the first line is the populations and the subsequent lines are the bi-allelic counts in each population for a number of SNPs. The first and second allele type has no meaning and can be chosen arbitrarily. 

## Running AdmixtureBayes

After installation, run AdmixtureBayes from an arbitrary folder by

```bash
$ AdmixtureBayes run --input_file allele_counts.txt --outgroup population4
```

This will calculate the covariance matrix saved in a new file called *covariance_and_multiplier.txt* and put results in a file called *result_mc3.csv*. To obtain the posterior distributions of the minimal topologies, topologies and number of admixtures run

```bash
$ AdmixtureBayes run --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt 
```

The result is a file called *posterior_distributions.csv*. From this file a couple of different plots can be constructed

```bash
$ AdmixtureBayes plot --plot consensus_trees --posterior_distribution_file posterior_distributions.csv
$ AdmixtureBayes plot --plot top_node_trees --posterior_distribution_file posterior_distributions.csv
$ AdmixtureBayes plot --plot top_trees --posterior_distribution_file posterior_distributions.csv
```

## Increasing number of populations

As more populations are added to the input file, the more steps it will take for the MCMC to converge. By default, the number of MCMC steps is 10,000 (= m\*n=50\*200) which is only suitable for datasets with 4 or fewer populations. To increase the number of steps to 1,000,000 which is often enough to analyze 10 populations, use the command  

```bash
$ AdmixtureBayes run --input_file big_allele_counts.txt --outgroup population10 --n 20000
```

It is also possible to stop the chain using a stopping criteria. The stopping criteria calculates the Effective Sample Size of different summaries and stops if all of them are above a certain threshold (default is 200). To use the on-the-fly MCMC stopping criteria, the program needs the [rwty package](https://cran.r-project.org/web/packages/rwty/index.html) in R.
```bash
$ R
...
> install.packages("rwty")
```
And run AdmixtureBayes with the command.

```bash
$ AdmixtureBayes run --input_file big_allele_counts.txt --outgroup population10 --n 500000 --stop_criteria
```


## If AdmixtureBayes could not be installed
In case AdmixtureBayes would not install, it can still be run as
```bash
$ python admixturebayes/__main__.py run --input_file my_data.txt
```
However, it will not use the accelerated c-implementation of the covariance matrix.