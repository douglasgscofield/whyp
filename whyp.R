# whyp.R
#
# Calculate parentage using mixed-assay likelihood (Smouse et al. 2012 Journal
# of Heredity).
#
# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Science Centre, Umeå, Sweden
#
# douglas.scofield@plantphys.umu.se
# douglasgscofield@gmail.com
#
# These statistical tools were developed in collaboration with Peter Smouse
# (Rutgers University) and Victoria Sork (UCLA) and were funded by U.S. National
# Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.
#
# Use as you see fit.  No warranty regarding this code is implied nor should be
# assumed.  Send bug reports etc. to one of the above email addresses.
#
.whypVersion = "0.01"
#
# CHANGELOG
#
# 0.01: First prerelease


z = try(source("readGenalex.R"), silent=TRUE)
if (class(z) == "try-error") {
  z = try(source("whyp_functions/readGenalex.R"), silent=TRUE)
  if (class(z) == "try-error") {
    stop("readGenalex.R must be available, download it at https://github.com/douglasgscofield/popgen")
  }
}
source("whyp_functions/whyp_read_data.R")
source("whyp_functions/whyp_alleles.R")
source("whyp_functions/whyp_likelihood.R")
source("whyp_functions/whyp_extract_likelihood.R")
source("whyp_functions/whyp_backward_probability.R")
source("whyp_functions/whyp_check.R")
source("whyp_functions/whyp_utility.R")

# Configuration variables

.max.matches.to.report <- 10

# missing allele indications in raw datafiles
missing.allele <- c("0", "", " ")
genalex.missing.allele <- "0"

# missing allele indication in indexed genotype matrices
missing.allele.index <- -1

# treat single missing allele as homozygote
missing.single.allele.avail.methods <- c("as.missing","as.homozygote","as.is")

#missing.single.allele.method <- "as.missing"
missing.single.allele.method <- match.arg("as.missing",
                      missing.single.allele.avail.methods)

# any loci to drop from the entire dataset?
.drop.loci <- NULL

# probability of a null genotyping error
.Pr.null <- 0.02

# use the seedling mismatch heuristic
.use.heur.mismatch <- TRUE
.heur.tissue.methods <- c("seedling","pericarp")
.report.heur.mismatch <- FALSE

# probability assigned to a heuristic mismatch
.Pr.heur.mismatch <- 0.01
.heur.methods <- c("if.no.match","any")

# data sources we can use to make a parental assignment
.data.sources <- c("pericarp.priority", "use.all", "pericarp", "seedling")
.poke.data.sources <- c("pericarp.priority", "pericarp", "seedling")
.poke.holes.assay.running <- FALSE
.joint.method <- "pericarp.priority"

# neg-log-likelihood window within the min value for which we report candidates
.L.window <- Inf

# window to use for diagnostics produced by extract.via.L
.included.candidate.window <- c( LOGFUNC(10^2) )

#.included.candidate.window <- c(LOGFUNC(exp(1)^3), 
#                LOGFUNC(10^1), 
#                LOGFUNC(10^2), 
#                LOGFUNC(10^3))

whyp <- function(dat, 
                 recruits=dat$recruit.names, 
                 heur.method="if.no.match", 
                 heur.tissue="seedling")
{
    df <- data.frame()
    i <- 0
    for (recruit in recruits) {
        LOD.recruit <- LOD.Xkl(Xkl=Xkl.recruit.single(recruit, dat),
                               dat=dat, 
                               heur.method=heur.method,
                               heur.tissue=heur.tissue)
        df <- rbind(df, LOD.package(LOD.scores=LOD.recruit))
        i <- i + 1
        if (i %% 20 == 0) 
            cat("done with recruit", i, recruit, "...\n")
    }
    cat("completed whyp analysis for", i, "recruits\n")
    ####
    df
}


extract.samples.with.n.missing.loci <- function(dat, n.missing.loci=NULL)
{
  # dat needs to have name-valued rows.  shouldn't be an issue
  if (! is.null(n.missing.loci)) {
    ml <- attr(dat,"n.missing.loci")
    ml.indexes <- which(ml == n.missing.loci)
    ans <- dat[ml.indexes,]
    ml <- ml[ml.indexes]
    for (att in names(attributes(dat))) 
      if (! att %in% c("dim","dimnames","names")) # copy all but these
        attr(ans,att) <- attr(dat,att)
    attr(ans,"n.missing.loci") <- ml  # update
    dat <- ans
  }
  dat
}


assemble.whyp.data <- function(mother=NULL, 
                              pericarp=NULL, 
                              seedling=NULL,
                              file.mother=NULL, 
                              file.pericarp=NULL, 
                              file.seedling=NULL,
                              report.foreign.recruit=NULL,
                              report.mismatch.recruit=NULL,
                              report.null.recruit=NULL,
                              file.allele.freqs=NULL, 
                              drop.loci=.drop.loci,
                              write.reports=TRUE)
{
  dat <- list()  # add it all to a big list
  is.mother.info = !is.null(file.mother) || !is.null(mother)
  is.seedling.info = !is.null(file.seedling) || !is.null(seedling)
  is.pericarp.info = !is.null(file.pericarp) || !is.null(pericarp)
  if (!is.mother.info) stop("must provide mother genotypes")
  if (is.pericarp.info && !is.seedling.info) {
    # for now, work around the inconvenience of having to provide both pericarp
    # and seedling genotypes... keep track in the list of which were actuall
    # provided
    file.seedling = file.pericarp
    seedling = pericarp
    attr(dat, "data.sources.given") <- c("pericarp")
  } else if (is.seeding.info && !is.pericarp.info) {
    file.pericarp = file.seedling
    pericarp = seedling
    attr(dat, "data.sources.given") <- c("seedling")
  } else if (!is.seeding.info && !is.pericarp.info) {
    stop("must provide pericarp or seedling genotypes, or both")
  } else {
    attr(dat, "data.sources.given") <- c("pericarp", "seedling")
  }
  mother <-   if (is.null(mother)) read.mother(file.mother) else mother
  pericarp <- if (is.null(pericarp)) read.pericarp(file.pericarp) else pericarp
  seedling <- if (is.null(seedling)) read.seedling(file.seedling) else seedling
  if (! is.null(drop.loci)) {
    cat("dropping loci",drop.loci," from mother, pericarp and seedling...\n")
    mother <- dropGenalexLoci(mother, drop.loci, quiet=TRUE)
    pericarp <- dropGenalexLoci(pericarp, drop.loci, quiet=TRUE)
    seedling <- dropGenalexLoci(seedling, drop.loci, quiet=TRUE)
  }
  if (TRUE) {
    # check consistency of loci used
    ml <- attr(mother, "locus.names"); nml <- length(ml)
    pl <- attr(pericarp, "locus.names"); npl <- length(pl)
    sl <- attr(seedling, "locus.names"); nsl <- length(sl)
    loci <- auto.drop.loci <- c()
    cat("# mother loci =", nml, "  # pericarp loci =", npl,
        "  # seedling loci =", nsl, "\n")
    if (nml == npl && all(ml == pl) && nml == nsl && all(ml == sl)) {
      cat("all mother, seedling, pericarp loci are identical and in same order\n")
      loci <- ml
      loci.identical <- TRUE
    } else {
      loci.identical <- FALSE
      if (all(pl %in% ml) && all(ml %in% pl)) {
        cat("mother and pericarp locus sets are identical, but not in same order\n")
        loci <- sort(ml)
      } else {
        cat("mother, pericarp, seedling locus sets are not identical\n")
        loci <- sort(intersect(ml, intersect(pl, sl)))
        auto.drop.loci <- setdiff(union(ml, union(pl, sl)), loci)
      }
    }
    if (length(auto.drop.loci) > 0) {
      cat("auto-dropping loci ", auto.drop.loci, " from all datasets\n")
      mother <- dropGenalexLoci(mother, auto.drop.loci, quiet=TRUE)
      pericarp <- dropGenalexLoci(pericarp, auto.drop.loci, quiet=TRUE)
      seedling <- dropGenalexLoci(seedling, auto.drop.loci, quiet=TRUE)
    }
    mother <- reorderGenalexLoci(mother, loci)
    pericarp <- reorderGenalexLoci(pericarp, loci)
    seedling <- reorderGenalexLoci(seedling, loci)
  }
  dat$mother <- mother
  cat("*** note: merge.recruit skipped to work around locus column reorder bug\n")
  # dat$recruit <- merge.recruit(pericarp=pericarp, seedling=seedling)
  dat$recruit <- list(pericarp=pericarp, seedling=seedling)
  if ("pericarp" %in% attr(dat, "data.sources.given")) {
      dat$recruit.names <- pericarp[,1]
  } else if ("seedling" %in% attr(dat, "data.sources.given")) {
      dat$recruit.names <- seedling[,1]
  } else stop("can't find recruit names")
  # do we read allele frequencies from a file, or calculate them from input data?
  if (! is.null(file.allele.freqs)) {
    dat$master.alist <- file.allele.freqs
  } else {
    # create allele frequencies from scratch
    cat("calculating allele frequencies from scratch...\n")
    mother.alist <- create.allele.list(dat$mother)
    pericarp.alist <- create.allele.list(dat$recruit$pericarp)
    seedling.alist <- create.allele.list(dat$recruit$seedling)
    recruit.alist <- merge.allele.list(pericarp.alist, seedling.alist)
    dat$foreign.alist <- find.foreign.alleles(mother.alist,
                                              foreign=recruit.alist,
                                              method="actual.counts")
    dat$master.alist <- merge.allele.list(mother.alist, dat$foreign.alist)
  }
  # create allele table and index genotypes using it
  cat("creating allele table and using it to index genotypes...\n")
  dat$allele.table <- create.allele.table(dat$master.alist)
  dat$imother <- index.genotypes(dat$mother, dat$allele.table)
  dat$ipericarp <- index.genotypes(dat$recruit$pericarp, dat$allele.table)
  dat$iseedling <- index.genotypes(dat$recruit$seedling, dat$allele.table)
  # accuracies of these depend on original and indexed versions sharing row order
  attr(dat$mother,"n.missing.loci") <- attr(dat$imother,"n.missing.loci")
  for (a in c("n.missing.loci","name.missing.loci")) {
    attr(dat$recruit$pericarp,a) <- attr(dat$ipericarp,a)
    attr(dat$recruit$seedling,a) <- attr(dat$iseedling,a)
  }
  if (! is.null(report.foreign.recruit) && ! is.null(dat$foreign.alist)
    && write.reports) {
    cat("writing report to file", report.foreign.recruit, "\n")
    sink(file=report.foreign.recruit)
    report.foreign.alleles(dat)
    sink()
  }
  # add names for check.recruit.mismatch()
  cat("checking for mismatches between recruit pericarp and seedling...\n")
  if (! is.null(report.mismatch.recruit)) {
    if (write.reports) {
      cat("writing report to file", report.mismatch.recruit, "\n")
      sink(file=report.mismatch.recruit)
      dat <- check.recruit.mismatch(dat, report=TRUE)
      sink()
    }
  } else {
    dat <- check.recruit.mismatch(dat, report=FALSE)
  }
  cat("editing mismatched genotypes...\n")
  dat <- edit.recruit.mismatch(dat, method="set.all.missing", update.attributes=TRUE)
  cat("checking for potential null genotyping errors between recruit pericarps and all mothers...\n")
  if (! is.null(report.null.recruit)) {
    if (write.reports) {
      cat("writing report to file", report.null.recruit, "\n")
      sink(file=report.null.recruit)
      dat <- check.recruit.null(dat, report=TRUE)
      sink()
    }
  } else {
    dat <- check.recruit.null(dat, report=FALSE)
  }
  #
  #
  dat
}


# old wrapper functions

full.mixed.assay <- function(read.freqs=FALSE)
{
  dat <- assemble.whyp.data(read.freqs=read.freqs)
  report.foreign.alleles(dat)
  Ljk.r <- mixed.assay(dat)
  print.mixed.assay(Ljk.r)
  invisible(Ljk.r)
}


mixed.assay <- function(dat, Pr.null=.Pr.null)
{
  # calculate recruit X mother X locus likelihoods
  cat("mixed.assay Pr.null =",Pr.null,"\n")
  Xjkl.p <- Xjkl.pericarp(dat$ipericarp, dat$imother, 
              dat$allele.table, Pr.null=Pr.null)
  Xjkl.s <- Xjkl.seedling(dat$iseedling, dat$imother, 
              dat$allele.table)
  # assemble recruit X mother likelihoods
  Ljk.r <- Ljk.recruit(Xjkl.pericarp=Xjkl.p, Xjkl.seedling=Xjkl.s)
  attr(Ljk.r, "data.sources.given") <- attr(dat, "data.sources.given")
  attr(Ljk.r, "data") <- dat
  ####
  ####
  Ljk.r
}

