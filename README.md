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

To make it you should be able to do:
```
./configure
make
```
For abbreviated information on the available options do
```
gsi_sim --help
```
and for compete information, try
```
gsi_sim --help-full
```



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

