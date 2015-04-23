# gsi_sim

The command line program `gsi_sim` is designed to do two separate tasks:
  1. Infer the population of origin and the mixing proportions of organisms
(typically fish like salmon) sampled from a mixed collection, whilst jointly
estimating the mixing proportions.  This is known in fisheries as 
"Genetic Stock Identification" or GSI
  2. Conduct simulations of such inference to predict the power of the genetic
markers assembled in a "genetic basline" to conduct.  

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
./configure
make
```

That creates the executable `gsisim`. Once you have done that. For abbreviated information on the available options do
```sh
./gsisim --help
```
and for compete information, try
```
./gsisim --help-full
```


## Using gsi_sim to find duplicate samples
One of the ways we use `gsi_sim` in our lab is to search for samples that have identical (or
near identical) genotypes.  Here is an example of how to do that on an example
data set included in the repository:
```
time ./gsisim -b test_data/geno-match-example-data.txt  --close-match-base 8 70 | grep CLOSE_MATCH_GENOTYPES

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

For more information on that do `./gsisim --help-full` and find the section on `--close-match-base`.


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

