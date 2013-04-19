Whyp (Who's Your Parent): Parentage Assignment
----------------------------------------------

Whyp provides R functions for maximum-likelihood assignment of parentage of
progeny based on genotypes from progeny and/or progeny-associated maternal
tissue (e.g. pericarps), and genotypes of potential parents.  Whyp handles two
kinds of genotyping error: null alleles (common in maternal tissue), and
mistyping (completely wrong genotype).

Development of whyp to date has focused on identifying the maternal parent
when potentially error-containing pericarp genotypes are available.

The whyp API is very rough at this stage.  Genotype input for parents,
pericarps and progeny is each in the form of an annotated `data.frame` produced by
`readGenalex()` (available at <https://github.com/douglasgscofield/popgen>).

  `assemble.whyp.data()` : load data for a whyp run into a large list

  `whyp()` : calculate parental assignments

The file `example.R` contains some example code which uses these functions.  The data files themselves are not included.  `assemble.v1()` and `assemble.v2()` load data for a whyp run.  `whyp.analysis()` runs a whyp analysis on a loaded dataset, and includes alternatives for specifying an analysis and returning results. 

The file `whyp.R` and the entire folder `whyp_functions/` are required.  See what I mean about a rough API? :-)

* * *

These statistical tools were developed in collaboration with Peter Smouse
(Rutgers University) and Victoria Sork (UCLA) and were funded by U.S. National
Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.

* * *

### Publications

Whyp (including earlier versions) has been used the following publications:

Smouse, P. E., V. L. Sork, D. G. Scofield, and D. Grivet. 2012. Using Seedling
and Pericarp Tissues to Determine Maternal Parentage of Dispersed Valley Oak
Recruits. _Journal of Heredity_ 103:250-259.

Scofield, D. G., V. L. Sork, and P. E. Smouse. 2010. Influence of acorn
woodpecker social behaviour on transport of coast live oak (<i>Quercus
agrifolia</i>) acorns in a southern California oak savanna. _Journal of Ecology_
98:561-572.

Scofield, D. G., V. R. Alfaro, V. L. Sork, D. Grivet, E. Martinez, J. Papp, A.
R. Pluess et al. 2011. Foraging patterns of acorn woodpeckers (<i>Melanerpes
formicivorus</i>) on valley oak (<i>Quercus lobata</i> Née) in two California oak
savanna-woodlands. _Oecologia_ 166:187-196.


