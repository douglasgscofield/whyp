# functions for calculating likelihoods

LOGFUNC.name <- "log10"
LOGFUNC <- function(x) { do.call(LOGFUNC.name,list(x)) }

# Pericarp priority: when pericarp is both present (not missing) and certain
# (the probability == 1 or log(probability) == 0), then make the seedling entry
# uninformative (the seedling probability = 1, which provides full support for
# loci that are certain in the pericarp)
#
# Use all: combined all pericarp and seedling information at once.  No
# distrimination among data sources on a locus-by-locus basis as for
# pericarp.priority.

Xkl.recruit.single <- function(recruit, dat, Pr.null=.Pr.null)
{
    Xkl.pericarp <- Xkl.pericarp.single(recruit, dat, Pr.null=Pr.null)
    # extend the pericarp missing-loci array along the mother axis
    Xkl.pericarp.missing <- array(FALSE, dim=dim(Xkl.pericarp))
    Xl.pericarp.missing <- attr(Xkl.pericarp,"Xl.missing")
    for (mom in 1:dim(Xkl.pericarp.missing)[2])
        Xkl.pericarp.missing[,mom,] <- Xl.pericarp.missing
    Xkl.seedling.priority <- Xkl.seedling <- Xkl.seedling.single(recruit, dat)
    # don't know why we verified that pericarp provided perfect information
    #Xkl.seedling.priority[Xkl.pericarp.missing == FALSE & Xkl.pericarp == 1] <- 1
    if (exists("retro.PuS") && retro.PuS) 
        Xkl.seedling.priority[Xkl.pericarp.missing == FALSE & Xkl.pericarp == 1] <- 1
    else 
        Xkl.seedling.priority[Xkl.pericarp.missing == FALSE] <- 1
    ans <- list(Xkl.pericarp=Xkl.pericarp,
                Xkl.seedling=Xkl.seedling,
                Xkl.seedling.priority=Xkl.seedling.priority,
                settings=list(Pr.null=Pr.null))
    class(ans) <- c("mixed.assay.probability.single", class(ans))
    ####
    ans
}


Lk.recruit.single <- function(recruit, dat, Pr.null=.Pr.null)
# Return likelihoods for the recruit by calling the Xkl routines
{
    Xkl <- Xkl.recruit.single(recruit=recruit, dat=dat, Pr.null=Pr.null)
    Lk.pericarp <- apply(-LOGFUNC(Xkl$Xkl.pericarp), c(1,2), sum)
    Lk.seedling <- apply(-LOGFUNC(Xkl$Xkl.seedling), c(1,2), sum)
    Lk.seedling.priority <- apply(-LOGFUNC(Xkl$Xkl.seedling.priority), c(1,2), sum)
    ans <- list(Lk.recruit.use.all=Lk.pericarp + Lk.seedling,
                Lk.recruit.priority=Lk.pericarp + Lk.seedling.priority,
                Lk.pericarp=Lk.pericarp, 
                Lk.seedling=Lk.seedling,
                Lk.seedling.priority=Lk.seedling.priority)
    class(ans) <- c("mixed.assay.single", class(ans))
    ans
}


#### #### #### #### #### #### #### ####


print.LOD <- function(LOD.scores, Lik.window=.included.candidate.window, sep="\t")
{
    df <- LOD.package(LOD.scores, Lik.window)
    print.data.frame(df)
}


LOD.package <- function(LOD.scores, Lik.window=.included.candidate.window)
{
    df <- data.frame(recruit="", n.loci=0, n.seed.loci=0, n.peri.loci=0,
                     settings="", method="", 
                     n.cand=0, 
                     name.cand="",
                     min.Lik.cand=0, 
                     n.min.Lik.cand=0, 
                     name.min.Lik.cand="", 
                     prob.random.cand=-1,
                     LOD.cand=0, 
                     delta.cand=0, 
                     delta.name="", 
                     heur.settings="",
                     heur.n.cand=-1,
                     heur.name.cand="",
                     heur.min.Lik.cand=-1,
                     heur.n.min.Lik.cand=-1,
                     heur.name.min.Lik.cand="",
                     heur.prob.random.cand=-1,
                     heur.LOD.cand=-1,
                     heur.delta.cand=-1,
                     heur.delta.name="",
                     window=0, 
                     n.window.cand=0,
                     heur.n.window.cand=0)[0,]
    settings <- unlist(LOD.scores$settings[names(LOD.scores$settings) != "logfunc"])
    fields.1 <- list(recruit=LOD.scores$recruit, n.loci=LOD.scores$n.loci$n.loci,
                     n.seed.loci=LOD.scores$n.loci$n.seed.loci,
                     n.peri.loci=LOD.scores$n.loci$n.peri.loci,
                     settings=paste(sep="=", collapse=";", names(settings), settings))
    for (m in names(LOD.scores$results)) {
        n.min.Lik.cand <- LOD.scores$results[[m]]$n.min.Lik.cand
        fields.2 <- list(method=m, 
                         n.cand=LOD.scores$results[[m]]$n.cand,
                         name.cand=LOD.scores$results[[m]]$name.cand,
                         min.Lik.cand=LOD.scores$results[[m]]$min.Lik.cand[1],
                         n.min.Lik.cand=n.min.Lik.cand,
                         name.min.Lik.cand=paste(collapse=";", 
                             names(LOD.scores$results[[m]]$min.Lik.cand)),
                         prob.random.cand=LOD.scores$results[[m]]$prob.random.recruit,
                         LOD.cand=LOD.scores$results[[m]]$LOD.cand[1],
                         # delta.cand is +1 to go past min cand(s) with itself
                         delta.cand=LOD.scores$results[[m]]$delta[n.min.Lik.cand+1],
                         delta.name=if (LOD.scores$results[[m]]$n.cand == 0)
                                        NA
                                    else
                                        names(LOD.scores$results[[m]]$delta[n.min.Lik.cand+1]),
                         heur.settings="",
                         heur.n.cand=-1,
                         heur.name.cand="",
                         heur.min.Lik.cand=-1,
                         heur.n.min.Lik.cand=-1,
                         heur.name.min.Lik.cand="",
                         heur.prob.random.cand=-1,
                         heur.LOD.cand=-1,
                         heur.delta.cand=-1,
                         heur.delta.name=""
                         )

        if (! is.null(LOD.scores$heur.results[[m]])) {
            heur.settings <- unlist(LOD.scores$heur.settings)
            fields.2$heur.settings=paste(sep="=", collapse=";", names(heur.settings), heur.settings)
            heur.n.min.Lik.cand <- LOD.scores$heur.results[[m]]$n.min.Lik.cand
            fields.2$heur.n.cand <- LOD.scores$heur.results[[m]]$n.cand
            fields.2$heur.name.cand <- LOD.scores$heur.results[[m]]$name.cand
            fields.2$heur.min.Lik.cand <- LOD.scores$heur.results[[m]]$min.Lik.cand[1]
            fields.2$heur.n.min.Lik.cand <- heur.n.min.Lik.cand
            fields.2$heur.name.min.Lik.cand <- paste(collapse=";", 
                                    names(LOD.scores$heur.results[[m]]$min.Lik.cand))
            fields.2$heur.prob.random.cand=LOD.scores$heur.results[[m]]$prob.random.recruit
            fields.2$heur.LOD.cand <- LOD.scores$heur.results[[m]]$LOD.cand[1]
            fields.2$heur.delta.cand <- LOD.scores$heur.results[[m]]$delta[heur.n.min.Lik.cand+1]
            fields.2$heur.delta.name <- if (LOD.scores$heur.results[[m]]$n.cand == 0)
                                            NA
                                        else
                                            names(LOD.scores$heur.results[[m]]$delta[heur.n.min.Lik.cand+1])
        }

        for (w in Lik.window) {
            .Lik.window.candidates <- function(w, n.cand="", results="") {
                if (fields.2[[n.cand]] == 0)
                    0
                else {
                    ans <- sum(LOD.scores[[results]][[m]]$delta <= w)
                    if (ans == length(LOD.scores[[results]][[m]]$delta)) {
                        if ("random" %in% unlist(strsplit(names(LOD.scores[[results]][[m]]$delta),":")))
                            ans <- ans - 1  # remove the random candidate
                        }
                        ans
                    }
            }
            fields.3 <- list(window=w, 
                             n.window.cand=.Lik.window.candidates(w, "n.cand", "results"),
                             heur.n.window.cand=if (is.null(LOD.scores$heur.results[[m]])) {
                                -1
                             } else { .Lik.window.candidates(w, "heur.n.cand", "heur.results") }
                            )
            ####
            new.assignment.row <- data.frame(c(fields.1, fields.2, fields.3))
            df <- rbind(df, new.assignment.row)
        }
    }
    rownames(df) <- 1:nrow(df)
    ####
    ####
    df
}


heur.update.Xkl <- function(Xkl.s, heur.mismatch, P, method=.heur.methods)
{
    if (.poke.holes.assay.running)
        return(NULL)  # for now
    method <- match.arg(method)
    if (method == "if.no.match" && any(P > 0))
        return(NULL)
    # list of index triples for the mismatch check
    mm <- cbind(1, 1:(dim(Xkl.s)[2]), heur.mismatch)
    # "disable" triples (via NA) where the mismatch locus is 0 (no actual mismatch)
    mm[mm[,3] == 0, 3] <- NA
    if (method == "if.no.match") {
        # "disable" triples where there's a match without mismatch
        mm[P > 0, 2] <- NA
    }
    # i've rethought this and it's better to know we can't find a match
    #if (all(is.na(apply(mm, 1, prod))) && method == "if.no.match")  # all triples include NA
    #    return(NULL)
    # apply the heuristic probability where there's a mismatch
    Xkl.s[mm] <- .Pr.heur.mismatch
    ####
    Xkl.s
}


LOD.Xkl <- function(Xkl, dat, method=.data.sources, heur.method="if.no.match", heur.tissue="seedling")
{
    method <- if (missing(method)) {
                  if (all(c("pericarp", "seedling") %in% attr(dat, "data.sources.given")))
                      c("pericarp.priority", "pericarp", "seedling")
                  else attr(dat, "data.sources.given")
              } else if (length(method) == 1)
                  match.arg(method)
              else method
    if (! "mixed.assay.probability.single" %in% class(Xkl)) stop("wrong class")
    heur.method <- match.arg(heur.method, .heur.methods)
    if (length(method) == 1)
        heur.tissue = method
    heur.tissue <- match.arg(heur.tissue, .heur.tissue.methods)

    Xkl <- Xkl.recruit.heur.mismatch(Xkl, dat)

    recruit <- dimnames(Xkl$Xkl.seedling)[[1]]
    recruit.row <- which(dat$recruit$seedling[,1] == recruit)
    mothers <- dimnames(Xkl$Xkl.seedling)[[2]]
    n.loci <- if ("pericarp" %in% method) 
        attr(dat$ipericarp, "n.loci")
    else attr(dat$iseedling, "n.loci")
    n.seed.loci <- if ("seedling" %in% method) 
            n.loci - attr(dat$iseedling, "n.missing.loci")[recruit.row]
        else 0
    n.peri.loci <- if ("pericarp" %in% method) 
            n.loci - attr(dat$ipericarp, "n.missing.loci")[recruit.row]
        else 0

    results <- heur.results <- list()

    for (m in method) {
        # must set prob.random.recruit, P, s.mismatch.1, P.mismatch.1
        P <- P.heur <- NULL

        if (m == "use.all") {
            prob.random.recruit <- prod(Pr.recruit.at.random(recruit, dat))
            P <- (apply(Xkl$Xkl.seedling, c(1,2), prod) * apply(Xkl$Xkl.pericarp, c(1,2), prod))
            if (heur.tissue == "seedling") {
                xkl <- heur.update.Xkl(Xkl$Xkl.seedling,
                                       attr(Xkl,"heur.mismatch.1.seedling"), P, 
                                       method=heur.method)
                if (! is.null(xkl))
                    P.heur <- (apply(xkl, c(1,2), prod) * apply(Xkl$Xkl.pericarp, c(1,2), prod))
            } else if (heur.tissue == "pericarp") {
                xkl <- heur.update.Xkl(Xkl$Xkl.pericarp,
                                       attr(Xkl,"heur.mismatch.1.pericarp"), P, 
                                       method=heur.method)
                if (! is.null(xkl))
                    P.heur <- (apply(xkl, c(1,2), prod) * apply(Xkl$Xkl.seedling, c(1,2), prod))
            }
        } else if (m == "pericarp.priority") {
            prob.random.recruit <- prod(Pr.recruit.at.random(recruit, dat))
            P <- (apply(Xkl$Xkl.seedling.priority, c(1,2), prod) * apply(Xkl$Xkl.pericarp, c(1,2), prod))
            if (heur.tissue == "seedling") {
                xkl <- heur.update.Xkl(Xkl$Xkl.seedling.priority,
                                       attr(Xkl,"heur.mismatch.1.seedling.priority"), P, 
                                       method=heur.method)
                if (! is.null(xkl))
                    P.heur <- (apply(xkl, c(1,2), prod) * apply(Xkl$Xkl.pericarp, c(1,2), prod))
            } else if (heur.tissue == "pericarp") {
                xkl <- heur.update.Xkl(Xkl$Xkl.pericarp,
                                       attr(Xkl,"heur.mismatch.1.pericarp.priority"), P, 
                                       method=heur.method)
                if (! is.null(xkl))
                    P.heur <- (apply(xkl, c(1,2), prod) * apply(Xkl$Xkl.seedling, c(1,2), prod))
            }
        } else if (m == "seedling") {
            prob.random.recruit <- prod(Pr.seedling.at.random(recruit, dat))
            P <- apply(Xkl$Xkl.seedling, c(1,2), prod)
            if (heur.tissue == "seedling") {
                xkl <- heur.update.Xkl(Xkl$Xkl.seedling,
                                       attr(Xkl,"heur.mismatch.1.seedling"), P, 
                                       method=heur.method)
                if (! is.null(xkl))
                    P.heur <- apply(xkl, c(1,2), prod)
            }
        } else if (m == "pericarp") {
            prob.random.recruit <- prod(Pr.pericarp.at.random(recruit, dat))
            P <- apply(Xkl$Xkl.pericarp, c(1,2), prod)
            if (heur.tissue == "pericarp") {
                xkl <- heur.update.Xkl(Xkl$Xkl.pericarp,
                                       attr(Xkl,"heur.mismatch.1.pericarp"), P, 
                                       method=heur.method)
                if (! is.null(xkl))
                    P.heur <- apply(xkl, c(1,2), prod)
            }
        } else stop("unrecognized method")
            
        calculate.candidates <- function(P) {
            Xk.cand <- rev(sort(P[P > 0]))  # sorted in decreasing Xk
            n.cand <- length(Xk.cand)
            #cat(recruit, "n candidates with method", m, "=", n.cand, "\n")
            if (n.cand == 0) {
                name.cand <- NA
                min.Lik.cand <- LOD.cand <- delta <- NA
                n.min.Lik.cand <- 0
            } else if (n.cand >= 1) {
                name.cand <- paste(collapse=";", names(Xk.cand))
                which.min.Lik <- which(Xk.cand == Xk.cand[1])
                min.Lik.cand <- -LOGFUNC(Xk.cand[which.min.Lik])
                n.min.Lik.cand <- length(which.min.Lik)
                LOD.cand <- LOGFUNC(Xk.cand[which.min.Lik] / prob.random.recruit)
                if (n.cand == 1) {
                    delta <- LOGFUNC(Xk.cand[1] / c(Xk.cand, prob.random.recruit))
                    names(delta) <- paste(sep=":", names(Xk.cand[1]), 
                                                c(names(Xk.cand[1]), "random"))
                } else if (n.cand >= 2) {
                    delta <- abs(LOGFUNC(Xk.cand[1]/Xk.cand))
                    names(delta) <- paste(sep=":",names(Xk.cand[1]),
                                                names(Xk.cand))
                    #delta <- abs(LOGFUNC(Xk.cand[1]/Xk.cand[-which.min.Lik]))
                    #names(delta) <- paste(sep=":",names(Xk.cand[1]),
                    #                              names(Xk.cand[-which.min.Lik]))
                }
            }
            list(Xk.cand=Xk.cand, n.cand=n.cand, name.cand=name.cand, min.Lik.cand=min.Lik.cand,
                 n.min.Lik.cand=n.min.Lik.cand, LOD.cand=LOD.cand, delta=delta)
        }

        names(P) <- mothers
        Cand <- calculate.candidates(P)
        results[[m]] <- list(method=m, n.cand=Cand$n.cand, name.cand=Cand$name.cand,
            min.Lik.cand=Cand$min.Lik.cand, n.min.Lik.cand=Cand$n.min.Lik.cand,
            LOD.cand=Cand$LOD.cand, delta=Cand$delta,
            prob.random.recruit=prob.random.recruit, Xk.cand=Cand$Xk.cand)

        if (! is.null(P.heur)) {
            names(P.heur) <- mothers
            Cand <- calculate.candidates(P.heur)
            heur.results[[m]] <- list(method=m, heur.method=heur.method, 
                n.cand=Cand$n.cand, name.cand=Cand$name.cand,
                min.Lik.cand=Cand$min.Lik.cand, n.min.Lik.cand=Cand$n.min.Lik.cand,
                LOD.cand=Cand$LOD.cand, delta=Cand$delta,
                prob.random.recruit=prob.random.recruit, Xk.cand=Cand$Xk.cand)
        }
    }

    ans <- list(recruit=recruit,
         settings=c(Xkl$settings, LOGFUNC.name=LOGFUNC.name, logfunc=LOGFUNC),
         n.loci=list(n.loci=n.loci, n.seed.loci=n.seed.loci, n.peri.loci=n.peri.loci),
         results=results, 
         heur.settings=c(heur.method=heur.method, heur.tissue=heur.tissue, heur.Pr=.Pr.heur.mismatch),
         heur.results=heur.results)
    class(ans) <- c("LOD", class(ans))
    ####
    ####
    ans
}


#### #### #### #### #### #### #### ####

Xkl.recruit.heur.mismatch <- function(Xkl, dat,
                                      report=.report.heur.mismatch,
                                      num.mismatch=1)
{

    # This heuristic looks across all mothers and checks for cases where there
    # is a match with the mother using the recruit pericarp, but there is
    # num.mismatch mismatched loci with the mother using the seedling (could be
    # seedling.priority subset).  If this case is found, and the heuristic is
    # enabled, then we create nonzero entry in the vector
    # heur.seedling.mismatch corresponding to the index of the seedling locus
    # at which the mismatch occurred.
   
    if (dim(Xkl$Xkl.pericarp)[1] != 1) stop("requires just one recruit at a time")
    xp <- Xkl$Xkl.pericarp
    xs <- Xkl$Xkl.seedling
    xsp <- Xkl$Xkl.seedling.priority
    recruit <- dimnames(xp)[[1]]
    mothers <- dimnames(xp)[[2]]
    heur.mismatch.seedling <- heur.mismatch.seedling.priority <- heur.mismatch.pericarp <- heur.mismatch.pericarp.priority <- rep(0, length(mothers))
    loci <- 1:dim(xp)[3]; 
    names(loci) <- dimnames(xp)[[3]]
    for (k in seq(along=mothers)) {
        if (prod(xp[1,k,loci]) > 0) {
            if (length(mm <- which(xs[1,k,loci] == 0)) == num.mismatch) {
                if (report) {
                    cat(sprintf("%s : mom %s (%d/%d) matches pericarp, mismatch seedling loc %s\n",
                                recruit, mothers[k], k, length(mothers), names(mm)))
                }
                heur.mismatch.seedling[k] <- mm
            }
            if (length(mm <- which(xsp[1,k,loci] == 0)) == num.mismatch) {
                if (report) {
                    cat(sprintf("%s : mom %s (%d/%d) matches pericarp, mismatch seedling.priority loc %s\n",
                                recruit, mothers[k], k, length(mothers), names(mm)))
                }
                heur.mismatch.seedling.priority[k] <- mm
            }
        }
        if (prod(xs[1,k,loci]) > 0) {
            if (length(mm <- which(xp[1,k,loci] == 0)) == num.mismatch) {
                if (report) {
                    cat(sprintf("%s : mom %s (%d/%d) matches seedling, mismatch pericarp loc %s\n",
                                recruit, mothers[k], k, length(mothers), names(mm)))
                }
                heur.mismatch.pericarp[k] <- mm
            }
        }
        if (prod(xsp[1,k,loci]) > 0) {
            if (length(mm <- which(xp[1,k,loci] == 0)) == num.mismatch) {
                if (report) {
                    cat(sprintf("%s : mom %s (%d/%d) matches seedling.priority, mismatch pericarp loc %s\n",
                                recruit, mothers[k], k, length(mothers), names(mm)))
                }
                heur.mismatch.pericarp.priority[k] <- mm
            }
        }
    }
    attr(Xkl, paste(sep=".","heur.mismatch",num.mismatch,"seedling")) <- heur.mismatch.seedling
    attr(Xkl, paste(sep=".","heur.mismatch",num.mismatch,"seedling.priority")) <- heur.mismatch.seedling.priority
    attr(Xkl, paste(sep=".","heur.mismatch",num.mismatch,"pericarp")) <- heur.mismatch.pericarp
    attr(Xkl, paste(sep=".","heur.mismatch",num.mismatch,"pericarp.priority")) <- heur.mismatch.pericarp.priority
    ####
    invisible(Xkl)
}


Ljk.recruit <- function(Xjkl.pericarp, Xjkl.seedling) 
{
    Ljkl.pericarp <- -log(Xjkl.pericarp)
    Xjl.pericarp.missing <- attr(Xjkl.pericarp,"Xjl.missing")
    # extend the pericarp missing-loci array along the mother axis
    Xjkl.pericarp.missing <- array(FALSE, dim=dim(Xjkl.pericarp))
    for (mom in 1:dim(Xjkl.pericarp.missing)[2])
        Xjkl.pericarp.missing[,mom,] <- Xjl.pericarp.missing
    Ljkl.seedling <- Ljkl.seedling.priority <- -log(Xjkl.seedling)
    # we currently do nothing special for missing seedlings
    Ljkl.seedling.priority[Xjkl.pericarp.missing == FALSE & Ljkl.pericarp == 0] <- 0
    Ljk.pericarp <- apply(Ljkl.pericarp, c(1,2), sum)
    Ljk.seedling <- apply(Ljkl.seedling, c(1,2), sum)
    Ljk.seedling.priority <- apply(Ljkl.seedling.priority, c(1,2), sum)
    ans <- list(Ljk.recruit.use.all=Ljk.pericarp + Ljk.seedling,
                Ljk.recruit.priority=Ljk.pericarp + Ljk.seedling.priority,
                Ljk.pericarp=Ljk.pericarp, 
                Ljk.seedling=Ljk.seedling,
                Ljk.seedling.priority=Ljk.seedling.priority)
    class(ans) <- c("mixed.assay", class(ans))
    ans
}


Xkl.seedling.single <- function(seedling, dat)
# produce a matrix of Xkl values for one seedling X mothers X loci
{
    if (!any("matrix" %in% class(dat$iseedling))) stop("iseedling not matrix")
    iseedling <- dat$iseedling[seedling,,drop=FALSE]
    imother <- dat$imother
    if (!any("matrix" %in% class(imother))) stop("imother not matrix")
    allele.table <- dat$allele.table
    if (!any("matrix" %in% class(dat$allele.table))) stop("allele.table not matrix")
    loci <- attr(dat$iseedling, "locus.names")
    n.loci <- attr(dat$iseedling, "n.loci")
    locus.columns <- attr(dat$iseedling, "locus.columns")
    mother.names <- rownames(imother); n.mother <- length(mother.names)

    # NOTE: this Xjl.missing code could be a source of errors, as it's a
    # missing-allele detection procedure that has nothing to do with earlier
    # missing-allele detection procedures.  It is based on the same data but
    # does not use the output from that.  This is because the poke-holes assay
    # code does not update the missing-allele attributes associated with the
    # genotype arrays, which would probably be inefficient.
    Xl.missing <- array(FALSE, dim=c(1, n.loci), dimnames=list(seedling, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        if (missing.allele.index == iseedling[1,l1] || 
            missing.allele.index == iseedling[1,l2])
            if ((missing.allele.index == iseedling[1,l1] &&
                 !missing.allele.index == iseedling[1,l2]) ||
                 (missing.allele.index == iseedling[1,l2] &&
                 !missing.allele.index == iseedling[1,l1])) {
                # only one is missing
                #cat("Xkl.seedling.single: for",seedling,"locus",locus,"only one is missing\n")
                if (missing.single.allele.method == "as.missing") {
                    Xl.missing[1,locus] <- TRUE
                } else if (missing.single.allele.method == "as.homozygote") {
                    if (iseedling[1,l1] == missing.allele.index)
                        iseedling[1,l1] = iseedling[1,l2]
                    else if (iseedling[1,l2] != missing.allele.index)
                        stop("missing single allele data inconsistency")
                    else
                        iseedling[1,l2] = iseedling[1,l1]
                } else stop("unknown missing single allele method")
            } else {
                Xl.missing[1,locus] <- TRUE
            }
    }
    Xkl <- array(0, dim=c(1, n.mother, n.loci), 
                 dimnames=list(seedling, mother.names, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        for (k in 1:n.mother) {
            m1 <- imother[k,l1]
            m2 <- imother[k,l2]
            mother.missing <- if (missing.allele.index %in% c(m1,m2)) TRUE 
                              else FALSE
            s1 <- iseedling[1,l1]
            s2 <- iseedling[1,l2]
            if (Xl.missing[1,locus]) {
                Prob <- 1  # could be anything
            } else if (mother.missing) {
                # Hardy-Weinberg proportion based on seedling alleles
                p <- (allele.table[loci[locus],s1] * 
                        allele.table[loci[locus],s2])
                Prob <- if (s1 != s2) 2*p 
                        else p
            } else {
                Prob <- 0
                if (s1 == s2) { # seedling homozygote
                    if (m1 == m2) { # mother homozygote
                        if (m1 != s1) { # homozygotes don't match
                            Prob <- 0  
                        } else { # homozygotes match, father gave matching allele
                            Prob <- allele.table[loci[locus],s1]
                        }
                    } else { # mother heterozygote
                        if (m1 == s1 || m2 == s1) { # one mother allele matches,
                            # father provided the other
                            Prob <- 0.5 * allele.table[loci[locus],s1]
                        } else { # no chance mother made seedling
                            Prob <- 0  
                        }
                    }
                } else { # seedling heterozygote
                    if (m1 == m2) { # mother homozygote
                        if (m1 == s1) { # father provided s2
                            Prob <- allele.table[loci[locus],s2]
                        } else if (m1 == s2) { # father provided s1
                            Prob <- allele.table[loci[locus],s1]
                        } else { # no chance mother made seedling
                            Prob <- 0  
                        }
                    } else { # mother heterozygote
                        if ((m1 == s1 && m2 == s2) ||
                            (m1 == s2 && m2 == s1)) { # matching heterozygotes,
                            # father provided the other allele, don't know which
                            Prob <- 0.5 * (allele.table[loci[locus],s1] + 
                                        allele.table[loci[locus],s2])
                        } else if (m1 == s1 || m2 == s1) { # father provided s2
                            Prob <- 0.5 * allele.table[loci[locus],s2]
                        } else if (m1 == s2 || m2 == s2) { # father provided s1
                            Prob <- 0.5 * allele.table[loci[locus],s1]
                        } else { # no chance mother made seedling
                            Prob <- 0  
                        }
                    }
                }
            }
            Xkl[1,k,locus] <- Prob
        }
    }
    attr(Xkl,"Xl.missing") <- Xl.missing
    Xkl
}


Xjkl.seedling <- function(iseedling, imother, allele.table)
# produce an array of Xjkl values for all seedlings X mothers X loci
{
    if (!any("matrix" %in% class(iseedling))) stop("iseedling not matrix")
    if (!any("matrix" %in% class(imother))) stop("imother not matrix")
    if (!any("matrix" %in% class(allele.table))) stop("allele.table not matrix")
    loci <- attr(iseedling, "locus.names"); n.loci <- attr(iseedling, "n.loci")
    locus.columns <- attr(iseedling, "locus.columns")
    seedling.names <- rownames(iseedling); n.seedling <- length(seedling.names)
    mother.names <- rownames(imother); n.mother <- length(mother.names)

    # NOTE: this Xjl.missing code could be a source of errors, as it's a
    # missing-allele detection procedure that has nothing to do with earlier
    # missing-allele detection procedures.  It is based on the same data but
    # does not use the output from that.  This is because the poke-holes assay
    # code does not update the missing-allele attributes associated with the
    # genotype arrays, which would probably be inefficient.
    Xjl.missing <- array(FALSE, dim=c(n.seedling, n.loci), 
                         dimnames=list(seedling.names, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        for (j in 1:n.seedling) {
            if (missing.allele.index == iseedling[j,l1] || 
                missing.allele.index == iseedling[j,l2])
                Xjl.missing[j,locus] <- TRUE
        }
    }
    Xjkl <- array(0, dim=c(n.seedling,n.mother,n.loci), 
                 dimnames=list(seedling.names, mother.names, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        for (k in 1:n.mother) {
            m1 <- imother[k,l1]
            m2 <- imother[k,l2]
            mother.missing <- if (missing.allele.index %in% c(m1,m2)) TRUE 
                              else FALSE
            for (j in 1:n.seedling) {
                s1 <- iseedling[j,l1]
                s2 <- iseedling[j,l2]
                if (Xjl.missing[j,locus]) {
                    Prob <- 1  # could be anything
                } else if (mother.missing) {
                    # Hardy-Weinberg proportion based on seedling alleles
                    p <- (allele.table[loci[locus],s1] * 
                          allele.table[loci[locus],s2])
                    Prob <- if (s1 != s2) 2*p 
                            else p
                } else {
                    Prob <- 0
                    if (s1 == s2) { # seedling homozygote
                        if (m1 == m2) { # mother homozygote
                            if (m1 != s1) { # homozygotes don't match
                                Prob <- 0  
                            } else { # homozygotes match, father gave matching allele
                                Prob <- allele.table[loci[locus],s1]
                            }
                        } else { # mother heterozygote
                            if (m1 == s1 || m2 == s1) { # one mother allele matches,
                                # father provided the other
                                Prob <- 0.5 * allele.table[loci[locus],s1]
                            } else { # no chance mother made seedling
                                Prob <- 0  
                            }
                        }
                    } else { # seedling heterozygote
                        if (m1 == m2) { # mother homozygote
                            if (m1 == s1) { # father provided s2
                                Prob <- allele.table[loci[locus],s2]
                            } else if (m1 == s2) { # father provided s1
                                Prob <- allele.table[loci[locus],s1]
                            } else { # no chance mother made seedling
                                Prob <- 0  
                            }
                        } else { # mother heterozygote
                            if ((m1 == s1 && m2 == s2) ||
                                (m1 == s2 && m2 == s1)) { # matching heterozygotes,
                                # father provided the other allele, don't know which
                                Prob <- 0.5 * (allele.table[loci[locus],s1] + 
                                            allele.table[loci[locus],s2])
                            } else if (m1 == s1 || m2 == s1) { # father provided s2
                                Prob <- 0.5 * allele.table[loci[locus],s2]
                            } else if (m1 == s2 || m2 == s2) { # father provided s1
                                Prob <- 0.5 * allele.table[loci[locus],s1]
                            } else { # no chance mother made seedling
                                Prob <- 0  
                            }
                        }
                    }
                }
                Xjkl[j,k,locus] <- Prob
            }
        }
    }
    attr(Xjkl,"Xjl.missing") <- Xjl.missing
    Xjkl
}


Xkl.pericarp.single <- function(pericarp, dat, Pr.null=.Pr.null)
# produce an array of Xkl values for one pericarp X mothers X loci
{
    if (!any("matrix" %in% class(dat$ipericarp))) stop("ipericarp not matrix")
    ipericarp <- dat$ipericarp[pericarp,,drop=FALSE]
    imother <- dat$imother
    if (!any("matrix" %in% class(imother))) stop("imother not matrix")
    allele.table <- dat$allele.table
    if (!any("matrix" %in% class(dat$allele.table))) stop("allele.table not matrix")
    loci <- attr(dat$iseedling, "locus.names")
    n.loci <- attr(dat$iseedling, "n.loci")
    locus.columns <- attr(dat$iseedling, "locus.columns")
    mother.names <- rownames(imother); n.mother <- length(mother.names)
    Pr.null <- rep(Pr.null, n.loci); names(Pr.null) <- loci

    # NOTE: this Xjl.missing code could be a source of errors, as it's a
    # missing-allele detection procedure that has nothing to do with earlier
    # missing-allele detection procedures.  It is based on the same data but
    # does not use the output from that.  This is because the poke-holes assay
    # code does not update the missing-allele attributes associated with the
    # genotype arrays, which would probably be inefficient.
    Xl.missing <- array(FALSE, dim=c(1, n.loci), dimnames=list(pericarp, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        if (missing.allele.index == ipericarp[1,l1] || 
            missing.allele.index == ipericarp[1,l2])
            if ((missing.allele.index == ipericarp[1,l1] &&
                 !missing.allele.index == ipericarp[1,l2]) ||
                 (missing.allele.index == ipericarp[1,l2] &&
                 !missing.allele.index == ipericarp[1,l1])) {
                # only one is missing
                #cat("Xkl.pericarp.single: for",pericarp,"locus",locus,"only one is missing\n")
                if (missing.single.allele.method == "as.missing") {
                    Xl.missing[1,locus] <- TRUE
                } else if (missing.single.allele.method == "as.homozygote") {
                    if (ipericarp[1,l1] == missing.allele.index)
                        ipericarp[1,l1] = ipericarp[1,l2]
                    else if (ipericarp[1,l2] != missing.allele.index)
                        stop("missing single allele data inconsistency")
                    else
                        ipericarp[1,l2] = ipericarp[1,l1]
                } else stop("unknown missing single allele method")
            } else {
                Xl.missing[1,locus] <- TRUE
            }
            #Xl.missing[1,locus] <- TRUE
    }
    Xkl <- array(0, dim=c(1, n.mother, n.loci), 
                 dimnames=list(pericarp, mother.names, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        for (k in 1:n.mother) {
            m1 <- imother[k,l1]
            m2 <- imother[k,l2]
            mother.missing <- if (missing.allele.index == m1 ||
                                  missing.allele.index == m2) TRUE 
                              else FALSE
            p1 <- ipericarp[1,l1]
            p2 <- ipericarp[1,l2]

            if (Xl.missing[1,locus]) {
                Prob <- 1  # could be anything
            } else if (mother.missing) {
                # Hardy-Weinberg proportion based on pericarp alleles
                p <- (allele.table[loci[locus],p1] * 
                        allele.table[loci[locus],p2])
                Prob <- if (p1 != p2) 2*p 
                        else p
            } else {
                Prob <- 0
                if (p1 == p2) { # pericarp homozygote
                    if (m1 == m2) { # mother homozygote
                        if (m1 != p1) { # homozygotes don't match
                            Prob <- 0  
                        } else { # homozygotes match, could be null in pericarp
                            Prob <- 1 - (Pr.null[locus] * 
                                        (1 - allele.table[loci[locus],p1]))
                        }
                    } else { # mother heterozygote
                        if (m1 == p1) { # one mother allele matches, could be null
                            # in pericarp
                            Prob <- Pr.null[locus] * allele.table[loci[locus],m2]
                        } else if (m2 == p1) { # other allele matches, could be null
                            # in pericarp
                            Prob <- Pr.null[locus] * allele.table[loci[locus],m1]
                        } else { # no chance mother made pericarp
                            Prob <- 0  
                        }
                    }
                } else { # pericarp heterozygote
                    if ((m1 == p1 && m2 == p2) ||
                        (m1 == p2 && m2 == p1)) { # matching heterozygotes
                        Prob <- 1
                    } else { # het mom doesn't match het pericarp
                        Prob <- 0  
                    }
                }
            }
            Xkl[1,k,locus] <- Prob
        }
    }
    attr(Xkl,"Xl.missing") <- Xl.missing
    Xkl
}


Xjkl.pericarp <- function(ipericarp, imother, allele.table, Pr.null=.Pr.null)
# produce an array of Xjkl values for all pericarps X mothers X loci
{
    if (!any("matrix" %in% class(ipericarp))) stop("ipericarp not matrix")
    if (!any("matrix" %in% class(imother))) stop("imother not matrix")
    if (!any("matrix" %in% class(allele.table))) stop("allele.table not matrix")
    loci <- attr(ipericarp, "locus.names"); n.loci <- attr(ipericarp, "n.loci")
    # dummy up a null probability for now
    Pr.null <- rep(Pr.null, n.loci); names(Pr.null) <- loci
    #Pr.null["Z11"] <- 0.05
    locus.columns <- attr(ipericarp, "locus.columns")
    pericarp.names <- rownames(ipericarp); n.pericarp <- length(pericarp.names)
    mother.names <- rownames(imother); n.mother <- length(mother.names)

    # NOTE: this Xjl.missing code could be a source of errors, as it's a
    # missing-allele detection procedure that has nothing to do with earlier
    # missing-allele detection procedures.  It is based on the same data but
    # does not use the output from that.  This is because the poke-holes assay
    # code does not update the missing-allele attributes associated with the
    # genotype arrays, which would probably be inefficient.
    Xjl.missing <- array(FALSE, dim=c(n.pericarp,n.loci), 
                         dimnames=list(pericarp.names, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        for (j in 1:n.pericarp) {
            if (missing.allele.index == ipericarp[j,l1] || 
                missing.allele.index == ipericarp[j,l2])
                Xjl.missing[j,locus] <- TRUE
        }
    }
    Xjkl <- array(0, dim=c(n.pericarp,n.mother,n.loci), 
                 dimnames=list(pericarp.names, mother.names, loci))
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]
        l2 <- l1 + 1
        for (k in 1:n.mother) {
            m1 <- imother[k,l1]
            m2 <- imother[k,l2]
            mother.missing <- if (missing.allele.index == m1 ||
                                  missing.allele.index == m2) TRUE 
                              else FALSE
            for (j in 1:n.pericarp) {
                p1 <- ipericarp[j,l1]
                p2 <- ipericarp[j,l2]

                if (Xjl.missing[j,locus]) {
                    Prob <- 1  # could be anything
                } else if (mother.missing) {
                    # Hardy-Weinberg proportion based on pericarp alleles
                    p <- (allele.table[loci[locus],p1] * 
                          allele.table[loci[locus],p2])
                    Prob <- if (p1 != p2) 2*p 
                            else p
                } else {
                    Prob <- 0
                    if (p1 == p2) { # pericarp homozygote
                        if (m1 == m2) { # mother homozygote
                            if (m1 != p1) { # homozygotes don't match
                                Prob <- 0  
                            } else { # homozygotes match, could be null in pericarp
                                Prob <- 1 - (Pr.null[locus] * 
                                            (1 - allele.table[loci[locus],p1]))
                            }
                        } else { # mother heterozygote
                            if (m1 == p1) { # one mother allele matches, could be null
                                # in pericarp
                                Prob <- Pr.null[locus] * allele.table[loci[locus],m2]
                            } else if (m2 == p1) { # other allele matches, could be null
                                # in pericarp
                                Prob <- Pr.null[locus] * allele.table[loci[locus],m1]
                            } else { # no chance mother made pericarp
                                Prob <- 0  
                            }
                        }
                    } else { # pericarp heterozygote
                        if ((m1 == p1 && m2 == p2) ||
                            (m1 == p2 && m2 == p1)) { # matching heterozygotes
                            Prob <- 1
                        } else { # het mom doesn't match het pericarp
                            Prob <- 0  
                        }
                    }
                }
                Xjkl[j,k,locus] <- Prob
            }
        }
    }
    attr(Xjkl,"Xjl.missing") <- Xjl.missing
    Xjkl
}


print.mixed.assay.single <- function(Lk, 
                                     data.source=.data.sources,
                                     include.min=TRUE, include.all=FALSE, 
                                     pretty=FALSE, print.Inf=FALSE, 
                                     diagnostics=FALSE,
                                     L.window=.L.window)
{
    if (abs(L.window) != L.window) stop("L.window must be >= 0 (may be +Inf)")
    recruit <- dimnames(Lk[[1]])[[1]]
    mothers <- dimnames(Lk[[1]])[[2]]
    cat("Recruit", recruit, "mixed-assay results across",
        length(mothers),"mothers\n\n")
    if (include.min) cat("min negLogL recruits, ")
    if (include.all) cat("all negLogL recruits, L.window =",L.window,", ")
    if (! print.Inf) cat("excluding Inf results, ")
    cat("\n\n")

    cat("recruit\tdata.source\tcriterion\tnegLogL\tcandidate\n")

    #cat(sep="","\n", recruit, " min negLogL\n")
    L.result <- list()
    L.result$pericarp.priority <- Lk[["Lk.recruit.priority"]][recruit,]
    L.result$use.all <- Lk[["Lk.recruit.use.all"]][recruit,]
    L.result$pericarp <- Lk[["Lk.pericarp"]][recruit,]
    L.result$seedling <- Lk[["Lk.seedling"]][recruit,]
    #L.seedling.priority <- Lk[["Lk.seedling.priority"]][recruit,]

    do.print <- function(row, name, data.source, L.window=L.window) {
        criterion <- if (L.window == 0) "min" else
                        if (L.window < Inf) sprintf("L.window=%G",L.window)
                        else "all"
        if (! print.Inf) {
            pos.Inf <- (row == +Inf)
            row <- row[! pos.Inf]
            if (sum(pos.Inf)) {
                row <- c(row, +Inf)
                names(row)[length(row)] <- paste("***",sum(pos.Inf),"candidates")
            }
        }
        min.L <- min(row)
        row <- row[row >= min.L & row <= (min.L + L.window)]
        row <- row[order(row, names(row))]
        for (i in seq(along=row)) {
            if (pretty)
                cat(sprintf("%-12s\t%-18s\t%s\t%10G\t%s\n", 
                            name, data.source, criterion, row[i], names(row)[i]))
            else
                cat(sprintf("%s\t%s\t%s\t%G\t%s\n", 
                            name, data.source, criterion, row[i], names(row)[i]))
        }
    }
    if (include.min) {
        for (ds in data.source)
            do.print(L.result[[ds]], recruit, ds, 0)
        if (pretty) cat("\n")
    }
    if (include.all) {
        for (ds in data.source)
            do.print(L.result[[ds]], recruit, ds, L.window=L.window)
        if (pretty) cat("\n")
    }
}


print.mixed.assay <- function(Ljk, 
                              data.source=.data.sources,
                              include.min=TRUE, include.all=FALSE, 
                              pretty=TRUE, print.Inf=FALSE, 
                              L.window=.L.window)
{
    if (! all(c("pericarp", "seedling") %in% attr(Ljk, "data.sources.given"))
        && missing(data.source)) {
        # data.sources.given was restricted somehow, and data.source argument is default
        data.source = attr(Ljk, "data.sources.given")
    }
    if (abs(L.window) != L.window) stop("L.window must be >= 0 (may be +Inf)")
    recruits <- dimnames(Ljk[[1]])[[1]]
    mothers <- dimnames(Ljk[[1]])[[2]]
    cat("Mixed-assay results\n\n")
    cat(length(recruits),"recruits, ",length(mothers),"mothers\n\n")
    if (include.min) cat("min negLogL recruits, ")
    if (include.all) cat("all negLogL recruits, L.window =",L.window,", ")
    if (! print.Inf) cat("excluding Inf results, ")
    cat("\n\n")

    cat("recruit\tdata.source\tcriterion\tnegLogL\tcandidate\n")

    for (recruit in recruits) {
        #cat(sep="","\n", recruit, " min negLogL\n")
        L.result <- list()
        L.result$pericarp.priority <- Ljk[["Ljk.recruit.priority"]][recruit,]
        L.result$use.all <- Ljk[["Ljk.recruit.use.all"]][recruit,]
        L.result$pericarp <- Ljk[["Ljk.pericarp"]][recruit,]
        L.result$seedling <- Ljk[["Ljk.seedling"]][recruit,]
        #L.seedling.priority <- Ljk[["Ljk.seedling.priority"]][recruit,]

        do.print <- function(row, name, data.source, L.window=L.window) {
            criterion <- if (L.window == 0) "min" else
                            if (L.window < Inf) sprintf("L.window=%G",L.window)
                            else "all"
            if (! print.Inf) {
                pos.Inf <- (row == +Inf)
                row <- row[! pos.Inf]
                if (sum(pos.Inf)) {
                    row <- c(row, +Inf)
                    names(row)[length(row)] <- paste("***",sum(pos.Inf),"candidates")
                }
            }
            min.L <- min(row)
            row <- row[row >= min.L & row <= (min.L + L.window)]
            row <- row[order(row, names(row))]
            for (i in seq(along=row)) {
                if (pretty)
                    cat(sprintf("%-12s\t%-18s\t%s\t%10G\t%s\n", 
                                name, data.source, criterion, row[i], names(row)[i]))
                else
                    cat(sprintf("%s\t%s\t%s\t%G\t%s\n", 
                                name, data.source, criterion, row[i], names(row)[i]))
            }
        }
        if (include.min) {
            for (ds in data.source)
                do.print(L.result[[ds]], recruit, ds, 0)
            if (pretty) cat("\n")
        }
        if (include.all) {
            for (ds in data.source)
                do.print(L.result[[ds]], recruit, ds, L.window=L.window)
            if (pretty) cat("\n")
        }
    }
}


annotate.Ljk <- function(Ljk.r, dat)
{
    recruits <- dimnames(Ljk.r$Ljk.pericarp)[[1]]
    mothers <- dimnames(Ljk.r$Ljk.pericarp)[[2]]
    attr(Ljk.r,"n.loci") <- attr(dat$ipericarp,"n.loci")
    pericarp.atts <- list(id=dimnames(dat$ipericarp)[[1]])
    if (! all(sort(pericarp.atts$id) == sort(recruits)))
        stop("pericarp and recruit ids do not match in Ljk and dat")
    seedling.atts <- list(id=dimnames(dat$iseedling)[[1]])
    if (! all(sort(pericarp.atts$id) == sort(recruits)))
        stop("seedling and recruit ids do not match in Ljk and dat")
    mother.atts <- list(id=dimnames(dat$imother)[[1]])
    if (! all(sort(mother.atts$id) == sort(mothers)))
        stop("mother ids do not match in Ljk and dat")
    for (att in c("n.missing.loci", "name.missing.loci",
                  "n.mismatch.loci", "name.mismatch.loci",
                  "n.double.missing.loci", "name.double.missing.loci")) {
        pericarp.atts[[att]] <- attr(dat$ipericarp,att)
        seedling.atts[[att]] <- attr(dat$iseedling,att)
    }
    for (att in c("n.missing.loci", "name.missing.loci")) {
        mother.atts[[att]] <- attr(dat$imother,att)
    }
    pericarp.atts <- as.data.frame(pericarp.atts,row.names=pericarp.atts$id)[recruits,]
    seedling.atts <- as.data.frame(seedling.atts,row.names=seedling.atts$id)[recruits,]
    mother.atts <- as.data.frame(mother.atts,row.names=mother.atts$id)[mothers,]
    for (att in c("locus.names","n.loci")) {
        attr(pericarp.atts,att) <- attr(dat$ipericarp,att)
        attr(seedling.atts,att) <- attr(dat$iseedling,att)
        attr(mother.atts,att) <- attr(dat$imother,att)
    }
    attr(Ljk.r,"annotations") <- list(pericarp=pericarp.atts, seedling=seedling.atts, mother=mother.atts)
    attr(Ljk.r,"data") <- dat
    #####
    Ljk.r
}


