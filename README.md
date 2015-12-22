# gsi_sim

The command line program `gsi_sim` is designed to do two separate tasks:
  1. Infer the population of origin and the mixing proportions of organisms
(typically fish like salmon) sampled from a mixed collection, whilst jointly
estimating the mixing proportions.  This is known in fisheries as 
"Genetic Stock Identification" or GSI
  2. Conduct simulations of such inference to predict the power of the genetic
markers assembled in a "genetic baseline".  

`gsi_sim` was originally written to do #2 above, but 
later was adopted to do  #1 as well.  

`gsi_sim` is written in C.  To get the source code from GitHub and compile it, do this
on the command line Unix/Linux or Mac machine (or Windows if you have set up a 
sane development environment on it):
```sh
# git the repo and submodules
git clone https://github.com/eriqande/gsi_sim.git
cd gsi_sim/
git submodule init
git submodule update

# make it
./Compile_gsi_sim.sh
```

That will compile up `gsi_sim-Darwin` on a Mac and `gsi_sim-Linux` on a Linux box.

Once you have done that. For abbreviated information on the available options do
```sh
./gsi_sim-Darwin --help
```
and for compete information, try
```
./gsi_sim-Darwin --help-full
```
replacing `Darwin` with `Linux` if you are on a Linux box and with `MINGW32_NT-6.1` if you are
on a Windows box.

## Executables

The executable files `gsi_sim-Darwin`, `gsi_sim-Linux`, `gsi_sim-MINGW32_NT-6.1` are provided as a courtesy, but are not guaranteed to have been compiled up from the latest commit. For that you should compile it up yourself (or see when the exectuable was last committed).
`gsi_sim-MINGW32_NT-6.1` is compiled by Eric using MINGW on a PC running as a virtual machine using VirtualBox on his Mac.  `gsi_sim-Linux` is compiled by Eric on our Ubuntu server at the lab.  `gsi_sim-Darwin` is compiled by Eric on his Mac laptop.

## Some Simple Examples
In the following examples, if you are running this on Linux, then you should replace
`gsi_sim-Darwin` with `gsi_sim-Linux`. If you are trying to run it on Windows then
replace `Darwin` with `MINGW32_NT-6.1` but be aware that you probably can't do the
processing of the output files via sed and awk on Windows, unless you have set that up.

### Self-assignment of baseline samples
Using a leave-one-out procedure, take all the baseline samples and assign them back to populations at the Unix command line.

```sh
# run the analysis and stick the output in a massive text file
 ./gsi_sim-Darwin -b test_data/snp349_baseline.txt --self-assign > dumpfile
 
# grab the resulting assignments out of that huge file using some unix/Linux tools
awk -F";" 'BEGIN {print "ID TopPop Score"} /SELF_ASSIGN_A_LA_GC_CSV:/ {print $1, $2, $3}' dumpfile | sed 's/SELF_ASSIGN_A_LA_GC_CSV:\///g;' > self-ass-results.txt 
```

After doing that we can read things into R (could have done that before, but R is not as fast as awk and sed for ripping
through large pure text files), and tally things up by population.

```r
library(dplyr)
library(stringr)

AssignmentTally <- read.table("self-ass-results.txt", header = TRUE) %>%
  tbl_df %>%
  mutate(FromPop = paste(str_split_fixed(ID, "_", 3)[,1], str_split_fixed(ID, "_", 3)[,2], sep ="_")) %>%
  group_by(FromPop, TopPop) %>%
  tally %>%
  arrange(FromPop, desc(n))
```

After that, `AssignmentTally` gives you many fish in the each group of the baseline were assigned to which populations. Here is
what it looks like when you do `as.data.frame(AssignmentTally)[1:25,]`

```
      FromPop           TopPop  n
1   F100_2014 BigQualicum_2014 83
2   F100_2014   Puntledge_2014 10
3   F100_2014        Inch_2006  1
4   F100_2014     Quinsam_2014  1
5   F103_2014    Capilano_2014 92
6   F103_2014     Quinsam_2014  2
7   F103_2014   Puntledge_2014  1
8   F104_2014   Robertson_2014 94
9   F105_2014   Puntledge_2014 85
10  F105_2014 BigQualicum_2014  8
11  F105_2014     Quinsam_2014  2
12  F106_2014     Quinsam_2014 92
13  F106_2014   Puntledge_2014  3
14  F132_2014      Salmon_2014 34
15 F1387_2014   Coldwater_2014 50
16  F150_2006        Inch_2006 63
17  F150_2006        Inch_2014  3
18  F150_2006    Chehalis_2014  2
19  F150_2006       Stave_2014  1
20  F150_2014        Inch_2014  2
21  F244_2014       Stave_2014 42
22  F244_2014    Chehalis_2014  2
23  F244_2014        Inch_2006  1
24  F244_2014        Inch_2014  1
25  M107_2014  Chilliwack_2014 93
```

### Doing MCMC to estimate mixture proportions
Here is another simple example.  Let's say we have the baseline `test_data/snp349_baseline.txt` and the
mixture file `test_data/snp349_baseline.txt`.  Then, to estimate the origin of each fish in the mixture
relative to the baseline using a Bayesian approach conditioned upon the baseline samples (e.g. not trying to
use the mixture to improve the baseline population allele frequencies), you can do like this:

```sh
 ./gsi_sim-Darwin -b test_data/snp349_baseline.txt  -t test_data/snp349_mixture.txt --mcmc-sweeps 25000 --mcmc-burnin 5000 > big_ol_output.txt 
```
It only takes a few seconds to do 30000 sweeps of MCMC.  Now, a few output files are produced
```
pop_pi_full_em_mle.txt           pop_pofz_one_step_em_mle.txt
pop_pi_posterior_means.txt       pop_pofz_posterior_means.txt
pop_pofz_full_em_mle.txt         pop_pofz_scaled_likelihoods.txt
```
"pofz" is the probability of the origin of each fish.  "pi" are the mixing proportions.  You get output both
from applying an EM algorithm to get an MLE and MCMC to get posterior means.  Crack those files open and check them
out.  They should be pretty self-explanatory.


### Compute Z-score to identify fish that may be from outside the baseline populations
Not yet written.

### Using gsi_sim to find duplicate samples
One of the ways we use `gsi_sim` in our lab is to search for samples that have identical (or
near identical) genotypes.  Here is an example of how to do that on an example
data set included in the repository (once again, replace `Darwin` with `Linux` if you are on a Linux box)
```
time ./gsi_sim-Darwin -b test_data/geno-match-example-data.txt  --close-match-base 8 70 | grep CLOSE_MATCH_GENOTYPES

# and the output for that is :
CLOSE_MATCH_GENOTYPES: T028958_?_10-6-2006 and T028959_?_10-6-2006 differ at 0 genotypes out of 95 non-missing loci
CLOSE_MATCH_GENOTYPES: T056350_M_9-29-2008 and T056351_M_9-29-2008 differ at 0 genotypes out of 85 non-missing loci

real	0m3.661s
user	0m3.639s
sys	0m0.019s

``` 
which shows that it found two pairs of samples that had identical genotypes out of about 4600 individuals,
and that required less than 4 seconds.

Output in this case also gets written to a file called `close_matches_baseline.txt` in the current 
working directory.

```sh
# look at the output:
cat close_matches_baseline.txt 

# Pairs of fish having 8 or fewer mismatching genotypes and at least 70 loci that are non-missing in each member of the pair 
# Each pair should be listed at most once. Indiv1 appears before Indiv2 in the data set.
Indiv1	Indiv2	NumMismatchLoci	NumNonMissingLoci
T028958_?_10-6-2006	T028959_?_10-6-2006	0	95
T056350_M_9-29-2008	T056351_M_9-29-2008	0	85
```

For more information on that do `./gsi_sim-Darwin --help-full` and find the section on `--close-match-base`.

The format of `geno-match-example-data.txt` should be pretty self evident.  It starts with the number of individuals
and the number of loci on the first line. Then a series of locus names, one on each line, and then a 
POP specifier with a population name after it.  Then a series of lines that each start with the individual
ID and are followed by the loci, with two columns for each diploid locus.  Alleles must be coded as integers.
We use A=1, C=2, G=3, and T=4 for SNPs.  Missing data at a locus is denoted with two zeroes (one for each allele).


## Citation

The main papers describing the simulation algorithms in gsi_sim are:
  1. Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. "An improved method for predicting the accuracy of genetic stock identification." Canadian Journal of Fisheries and Aquatic Sciences 65, no. 7 (2008): 1475-1486.
  2. Anderson, E. C. "Assessing the power of informative subsets of loci for population assignment: standard methods are upwardly biased." Molecular ecology resources 10, no. 4 (2010): 701-710.

There is no paper directly addressing the inference portion of gsi_sim, but it has seen use and citation in, for example:
  1. Israel, Joshua A., K. Jun Bando, Eric C. Anderson, and Bernie May. "Polyploid microsatellite data reveal stock complexity among estuarine North American green sturgeon (Acipenser medirostris)." Canadian Journal of Fisheries and Aquatic Sciences 66, no. 9 (2009): 1491-1504.
  2. Satterthwaite, William H., Michael S. Mohr, Michael R. O’Farrell, Eric C. Anderson, Michael A. Banks, Sarah J. Bates, M. Renee Bellinger et al. "Use of Genetic Stock Identification Data for Comparison of the Ocean Spatial Distribution, Size at Age, and Fishery Exposure of an Untagged Stock and Its Indicator: California Coastal versus Klamath River Chinook Salmon." Transactions of the American Fisheries Society 143, no. 1 (2014): 117-133.


## Terms 

As a work of the United States Government, this package is in the
public domain within the United States. Additionally, we waive
copyright and related rights in the work worldwide through the CC0 1.0
Universal public domain dedication.

See TERMS.md for more information.

