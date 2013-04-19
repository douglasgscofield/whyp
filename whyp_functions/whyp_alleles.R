# allele-related functions

report.foreign.alleles <- function(dat)
{
  if (is.null(dat$foreign.alist)) {
    warning("cannot report foreign alleles, no foreign allele list available")
    return
  }
  # report pericarps and progeny having non-mother alleles
  report.allele.carriers(dat$recruit$pericarp, dat$foreign.alist, 
                         "foreign pericarp")
  report.allele.carriers(dat$recruit$seedling, dat$foreign.alist, 
                         "foreign seedling")
}


report.allele.carriers <- function(dat, alist, note="")
{
  for (locus in names(alist)) {
    locus.data <- getGenalexLocus(dat, locus)
    alleles <- names(alist[[locus]])
    rows <- which(apply(locus.data, 1, function(x) any(alleles %in% x)))
    if (length(rows)) {
      cat(note, "allele carriers, locus",locus,"allele(s)",paste(alleles),"\n")
      printGenalexGenotype(dat, rows, callout.locus=locus)
    }
    cat("\n")
  }
}


report.allele.stats <- function(alist, round.to=c(4))
{
  funcs <- c("effective.n.alleles", "probability.identity", "sum")
  for (f in funcs) {
    for (r in round.to) {
      cat(f, "rounded to", r, "digits\n")
      print(round(unlist(lapply(alist, 
                                function(x)do.call(f, list(unlist(x))))), r))
    }
  }
}


# calculate allele counts for each locus from data 
create.allele.list <- function(data)
{
  loci <- attr(data,"locus.names")
  allele.list <- list() # missing.list <- list()
  for (locus in loci) {
    alleles <- as.vector(as.matrix(getGenalexLocus(data, locus)))
    missing <- alleles %in% missing.allele
    alleles <- alleles[! missing]
    tab <- table(alleles)
    allele.list[[locus]] <- as.list(tab)
    # missing.list[[locus]] <- sum(missing)
  }
  # if (max(unlist(missing.list)) > 0)
  #   attr(allele.list,"missing.list") <- missing.list
  allele.list
}


find.foreign.alleles <- function(alist, 
                                 foreign, 
                                 method=c("as.1", "actual.counts"))
{
  method <- match.arg(method)
  if (missing(foreign)) 
    stop("must supply foreign as source of novel alleles")
  result.alist <- list()
  for (locus in names(foreign)) {
    if (! locus %in% names(alist)) 
      stop("locus in foreign not in master")
    foreign.names <- names(foreign[[locus]])
    foreign.alleles <- foreign[[locus]][which(! foreign.names %in% 
                                                names(alist[[locus]]))]
    # cat(locus, "foreign.alleles:", paste(names(foreign.alleles)), "\n")
    if (length(foreign.alleles)) {  # alleles in foreign that are not in alist
      if (method == "as.1")
        result.alist[[locus]] <- lapply(foreign.alleles, function(x) 1)
      else if (method == "actual.counts")
        result.alist[[locus]] <- foreign.alleles
    }
  }
  ####
  result.alist
}


merge.allele.list <- function(alist1, alist2)
{
  # treat alist1 as our standard, add alist2 to it
  for (locus in names(alist2)) {
    if (! locus %in% names(alist1)) 
      alist1[[locus]] <- list()
    for (allele in names(alist2[[locus]])) {
      if (! allele %in% names(alist1[[locus]])) 
        alist1[[locus]][[allele]] <- alist2[[locus]][[allele]]
      else alist1[[locus]][[allele]] <- alist1[[locus]][[allele]] +
                                        alist2[[locus]][[allele]]
    }
    alist1[[locus]] <- alist1[[locus]][order(names(alist1[[locus]]))]
  }
  alist1
}


create.allele.table <- function(alist)
{
  loci <- sort(names(alist))
  max.alleles <- max(sapply(alist,length))
  # now turn this into a table of frequencies indexed by locus and allele
  # numbers
  allele.table <- matrix(0, nrow=length(loci), ncol=max.alleles)
  name.table <- aname.table <- matrix("", nrow=length(loci), ncol=max.alleles)
  rownames(allele.table) <- loci
  rownames(name.table) <- rownames(aname.table) <- loci
  frequency.list <- list()
  tmp <- 1:length(loci)
  names(tmp) <- loci
  locus.index.list <- as.list(tmp)
  index.list <- list()
  for (locus in loci) {
    alleles <- unlist(alist[[locus]])
    alleles <- alleles / sum(alleles)
    n.alleles <- length(alleles)
    allele.table[locus, 1:n.alleles] <- alleles
    name.table[locus, 1:n.alleles] <- paste(sep="/", locus, names(alleles))
    aname.table[locus, 1:n.alleles] <- names(alleles)
    frequency.list[[locus]] <- as.list(alleles)
    tmp <- 1:n.alleles
    names(tmp) <- names(alleles)
    index.list[[locus]] <- as.list(tmp)
  }
  attr(allele.table,"locus.names") <- loci
  attr(allele.table,"name.table") <- name.table
  attr(allele.table,"aname.table") <- aname.table
  attr(allele.table,"frequency.list") <- frequency.list
  attr(allele.table,"locus.index.list") <- locus.index.list
  attr(allele.table,"index.list") <- index.list
  ####
  allele.table
}


write.allele.freqs <- function(alist, file="allelefreqs.txt")
{
  write.table(list.to.data.frame(alist), file=file, sep="\t", quote=FALSE, 
              row.names=FALSE)
}


read.allele.freqs <- function(file="allelefreqs.txt")
{
  # locus<tab>allele<tab>freq
  data <- read.delim(file)
  if (any(names(data) != c("locus", "allele", "freq"))) 
    stop("headers wrong")
  # now make this an allele list 
  ####
  data.frame.to.list(data, by=c("locus", "allele"))
}


# calculate allele counts for each locus from data 
create.allele.list <- function(data)
{
  loci <- attr(data,"locus.names")
  allele.list <- list() # missing.list <- list()
  for (locus in loci) {
    alleles <- as.vector(as.matrix(getGenalexLocus(data, locus)))
    missing <- alleles %in% missing.allele
    alleles <- alleles[! missing]
    tab <- table(alleles)
    allele.list[[locus]] <- as.list(tab)
    # missing.list[[locus]] <- sum(missing)
  }
  attr(allele.list,"locus.names") <- loci
  # attr(allele.list,"missing.list") <- missing.list
  ####
  allele.list
}


effective.n.alleles <- function(x, count=TRUE) 
{
  if (count) 
    x <- x / sum(x)
  ####
  1 / sum(x*x)
}


pi.drop <- function(alist)
{
  pi.locus <- unlist(lapply(alist, function(x) do.call("probability.identity",
                                                       list(unlist(x)))))
  n.loci <- length(pi.locus)
  pi <- prod(pi.locus)
  g <- matrix(0, n.loci, n.loci-1)
  for (l1 in 1:n.loci) {
    g[l1,] <- pi.locus[-l1]
  }
  cat("mean dropping 1 =", mean(apply(g, 1, prod)), "\n")
  cat("stderr dropping 1 =", stderr(apply(g, 1, prod), na.rm=F), "\n")
  cat("n dropping 1 =", nrow(g), "\n")
  gg <- matrix(0, (n.loci * (n.loci - 1)), (n.loci - 2))
  for (l1 in 1:n.loci) {
    for (l2 in 1:(n.loci - 1)) {
      r.gg <- ((l1 - 1) * (n.loci - 1)) + l2
      cat("r.gg =", r.gg, " col g to drop =", -l2, "\n")
      gg[r.gg, ] <- g[l1, -l2]
    }
  }
  cat("mean dropping 2 =", mean(apply(gg, 1, prod)), "\n")
  cat("stderr dropping 2 =", stderr(apply(gg, 1, prod), na.rm=F), "\n")
  cat("n dropping 2 =", nrow(gg), "\n")
  ####
  return(gg)
}


probability.identity <- function(x, count=TRUE) 
{
  if (count) 
    x <- x/sum(x)

  # below, p_{i} is the frequency of the i-th allele at the locus

  # Formula from GenAlEx 6.1 Appendix 1, *** WHICH HAD TYPO ***
  #prob.identity.genalex <- (sum(x*x))^2 - sum(x*x*x*x)
  #cat("Probability Identity GenAlEx6 App. 1 =", prob.identity.genalex, "\n")

  # formula from Hedrick, Genetics of Populations, 4th ed.
  # PI = sum_{i} (p_{i}^4) + sum_{i} sum_{i < j} ((2*p_{i}*p_{j})^2)
  # which is probability of identity for homozygote draws plus probability
  # of identity for heterozygote draws
  het <- 0
  n <- length(x)
  if (n >= 2) {
    for (i in 1:(n-1))
      for (j in (i+1):n)
        het <- het + ((2 * x[i] * x[j]) ^ 2)
  }
  prob.identity.hedrick <- het + sum(x * x * x * x)
  # cat("Prob Identity Hedrick Eq. (11.12a) =", prob.identity.hedrick, "\n")
  ####
  prob.identity.hedrick
}


