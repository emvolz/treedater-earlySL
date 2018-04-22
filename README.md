# treedater-earlySL

This is a reanalysis of Ebola virus sequences originally presented by [Gire et al.](http://dx.doi.org/10.1126/science.1259657) and later reanalysied by [Moller et al.](https://doi.org/10.1073/pnas.1713314115).
Moller et al. suggest that substitution rate estimates based on these data are very sensitive to the population genetic prior used in Bayesian analyses. 
We're therefore curious what we estimate when we use a pseudo-maximum likelihood approach that does not incorporate such a prior.

## treedater 
The _treedater_ algorithm is [described here](https://doi.org/10.1093/ve/vex025). This is an [R package](https://github.com/emvolz/treedater), and we can very quickly estimate rates and dates from the 81 whole genomes used by Gire and Moller et al. 

Load the packages and data: 

```r
require(treedater)
require(ape)

( aln <- read.dna( 'gire_coding.fas', 'fasta') )
```

```
## 81 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 14517 
## 
## Labels:
## EBOV|KM034563|G3687|SierraLeone_G|2014-05-28
## EBOV|KM034562|G3686|SierraLeone_G|2014-05-28
## EBOV|KM034561|G3683|SierraLeone_G|2014-05-28
## EBOV|KM034560|G3682|SierraLeone_G|2014-05-28
## EBOV|KM034559|G3680|SierraLeone_G|2014-05-28
## EBOV|KM034558|G3679|SierraLeone_G|2014-05-28
## ...
## 
## Base composition:
##     a     c     g     t 
## 0.308 0.224 0.211 0.257
```

We will do a *very lazy* phylogenetic analysis. Afterall, it's Sunday and we need to watch _Lost in Space_. 
Following Moller, we will use a HKY model with no rate variation between sites. Then we will estimate a NJ tree. 

```r
D <- dist.dna( aln, model = 'F84', pairwise.del = TRUE )
tre <- unroot( bionj( D ) )
```


To estimate a molecular clock, we must extract the date from the sequence labels. These take the form of (for example) "EBOV|AA000000|EM104|SierraLeone_EM|2014-06-02"

```r
require(magrittr)
# note: times will be in units of years after Jan 1, 2014
sts <- strsplit( tre$tip.label, '\\|') %>% sapply( function(x) tail(x,1) ) %>% as.Date %>% -as.Date('2014-01-01' )/365
# note: making this a named vector corresponding to tree tip labels 
sts <- sts %>% as.numeric %>% setNames( tre$tip.label )
```

Now we can run _treedater_ with a relaxed clock

```r
( dtr.relaxed <- dater( tre, sts , s = 14e3 ) )
```



```
## 
## Phylogenetic tree with 81 tips and 80 internal nodes.
## 
## Tip labels:
## 	EBOV|KM034563|G3687|SierraLeone_G|2014-05-28, EBOV|KM034562|G3686|SierraLeone_G|2014-05-28, EBOV|KM034561|G3683|SierraLeone_G|2014-05-28, EBOV|KM034560|G3682|SierraLeone_G|2014-05-28, EBOV|KM034559|G3680|SierraLeone_G|2014-05-28, EBOV|KM034558|G3679|SierraLeone_G|2014-05-28, ...
## 
## Rooted; includes branch lengths.
## 
##  Time of common ancestor 
## 0.0906783357248483 
## 
##  Time to common ancestor (before most recent sample) 
## 0.369595636877891 
## 
##  Mean substitution rate 
## 0.000797488070378949 
## 
##  Strict or relaxed clock 
## relaxed 
## 
##  Coefficient of variation of rates 
## 0.279517382692339
```
And parametric bootstrap for the confidence intervals

```r
(boot.dtr.relaxed <-  parboot.treedater( dtr.relaxed , ncpu = 8)  )
```



```
##                            pseudo ML         2.5 %     97.5 %
## Time of common ancestor 0.0906783357 -2.726025e+01 0.18877756
## Mean substitution rate  0.0007974881  5.523629e-05 0.01151394
## 
##  For more detailed output, $trees provides a list of each fit to each simulation
```

And likewise with a strict clock, which is more like what was done in Moller et al

```r
( dtr.strict <- dater( tre, sts , s = 14e3 , strict=TRUE) )
```



```
## 
## Phylogenetic tree with 80 tips and 79 internal nodes.
## 
## Tip labels:
## 	EBOV|KM034563|G3687|SierraLeone_G|2014-05-28, EBOV|KM034562|G3686|SierraLeone_G|2014-05-28, EBOV|KM034561|G3683|SierraLeone_G|2014-05-28, EBOV|KM034560|G3682|SierraLeone_G|2014-05-28, EBOV|KM034559|G3680|SierraLeone_G|2014-05-28, EBOV|KM034558|G3679|SierraLeone_G|2014-05-28, ...
## 
## Rooted; includes branch lengths.
## 
##  Time of common ancestor 
## 0.0671399234498273 
## 
##  Time to common ancestor (before most recent sample) 
## 0.393134049152912 
## 
##  Mean substitution rate 
## 0.00101889360162669 
## 
##  Strict or relaxed clock 
## strict 
## 
##  Coefficient of variation of rates 
## 0
```

```r
(boot.dtr.strict <-  parboot.treedater( dtr.strict , ncpu = 8)  )
```
```
##                           pseudo ML         2.5 %      97.5 %
## Time of common ancestor 0.067139923 -3.7643027622 0.174085965
## Mean substitution rate  0.001018894  0.0001809049 0.005738618
## 
##  For more detailed output, $trees provides a list of each fit to each simulation
```

## Summary 
* Our rates are a little lower, but close to, what was reported by Moller et al. using the structured coalescent prior, especially with a strict clock.
* The estimated TMRCA is in early 2014, which is an acceptable estimate. This is intermediate between the estimates using structured coalescent and other priors reported by Moller. 
* We can't very well estimate the lower bound of TMRCA from this small aligment using the parametric bootstrap. The priors used in the Bayesian methods probably have a large influence on that. We also have overly wide confidence intervals. This could be improved by using a more time consuming bootstrap procedure. 
* Estimates could certainly be improved by doing a proper phylogenetic analysis and fine-tuning parameters
* All of the preceding analysis took less than one minute of computation time. 
