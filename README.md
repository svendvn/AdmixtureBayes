# AdmixtureBayes
## Purpose
AdmixtureBayes estimates posterior samples of admixture graphs (phylogenies incorporating admixture events) given an allele count file.

AdmixtureBayes is also shipped with commands to summarize and plot the results

## Installation


### Linux 64-bit

For Linux 64-bit, install the program in a python 2.7 environment using the commands
```bash
$ git clone https://github.com/svendvn/AdmixtureBayes
$ pip install dist/AdmixtureBayes-0.3-cp27-cp27mu-linux_x86_64.whl
```
It will install the necessary python dependencies. To se the plotting tools, the program dot (from grapvhiz) has to be accessible from the command line. Even if it is installed as one of the depencies it may not be visible to the command line. If so, install it manually
```bash
$ apt-get install graphviz
```


### Other OS

For other than linux 64-bit, AdmixtureBayes can be install by running the commands
```bash
$ git clone https://github.com/svendvn/AdmixtureBayes
$ python setup.py install
```

in a python 2.7 environment which will compile a C-file. Alternatively, AdmixtureBayes can be installed without the C-accelerated likelihood in a python 2.7 environment

```bash
$ git clone https://github.com/svendvn/AdmixtureBayes
$ pip install dist/AdmixtureBayes-0.3-py2-none-any.whl
```

It may also be necessary to install graphviz separately.

### Test installation

A test script is found in the *example/* folder together with a test dataset.
 
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

After installation, run AdmixtureBayes from any folder (containing the input file, here *allele_counts.txt* which is included in the folder *example/*).

```bash
$ AdmixtureBayes run --input_file big_allele_counts.txt --outgroup out
```

This will calculate the covariance matrix saved in a new file called *covariance_and_multiplier.txt* and put results in a file called *result_mc3.csv*. To obtain the posterior distributions of the minimal topologies, topologies and number of admixtures run

```bash
$ AdmixtureBayes posterior --input_file result_mc3.csv --covariance_matrix_file covariance_and_multiplier.txt 
```

The result is a file called *posterior_distributions.csv*. From this file a couple of different plots can be constructed

```bash
$ AdmixtureBayes plot --plot consensus_trees --posterior_distribution_file posterior_distributions.csv
$ AdmixtureBayes plot --plot top_node_trees --posterior_distribution_file posterior_distributions.csv
$ AdmixtureBayes plot --plot top_trees --posterior_distribution_file posterior_distributions.csv
$ AdmixtureBayes plot --plot estimates --posterior_distribution_file posterior_distributions.csv
```

## Increasing number of populations

As more populations are added to the input file, the more steps it will take for the MCMC to converge. By default there are 50 MCMC steps between each MCMCMC flip and the number of MCMCMC flips is 200. The total number of MCMC steps is therefore 200\*50=10,000 which is only suitable for datasets with 4 or fewer populations. To increase the number of steps to 1,000,000 which is often enough to analyze 10 populations, use the command  

```bash
$ AdmixtureBayes run --input_file big_allele_counts.txt --outgroup population10 --n 20000
```

It is also possible to stop the chain using a stopping criteria. The stopping criteria calculates the Effective Sample Size of different summaries and stops if all of them are above a certain threshold (default is 200). To do so, AdmixtureBayes calls the [rwty package](https://cran.r-project.org/web/packages/rwty/index.html) in R with the command Rscript. To use the stopping criteria it is therefore necessary to install rwty
```bash
$ R
...
> install.packages("rwty")
```
and run AdmixtureBayes with the command

```bash
$ AdmixtureBayes run --input_file big_allele_counts.txt --outgroup population10 --n 50000 --stop_criteria
```

The parameter --n then defines an upper limit on the number of iterations. 


## If AdmixtureBayes could not be installed
In case AdmixtureBayes would not install, it can still be run as
```bash
$ python admixturebayes/__main__.py run --input_file my_data.txt
```
However, it will not use the accelerated c-implementation of the covariance matrix.
