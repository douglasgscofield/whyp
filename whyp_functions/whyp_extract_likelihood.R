# functions for extracting according to likelihood criteria

extract.via.L <- function(Ljk.r, data.source=.data.sources,
                                 method=function(x) (x == min(x)),
                                 return.full.match=FALSE,
                                 drop.Inf=TRUE)
{
    # Returns a data frame with recruits extracted from Ljk.r according to their
    # negLogLik, as specified by the method function.  If there is more than one
    # mother, the mothers are separated by comma.

    data.source <- match.arg(data.source)
    recruits <- dimnames(Ljk.r[[1]])[[1]]
    mothers <- dimnames(Ljk.r[[1]])[[2]]
    label <- switch(data.source, pericarp.priority = "Ljk.recruit.priority",
                                 use.all = "Ljk.recruit.use.all",
                                 pericarp = "Ljk.pericarp",
                                 seedling = "Ljk.seedling")
    ans <- data.frame(recruit="", n.mother=0, mother="", data.source="", negLogL=0, 
                      candidate.window="", n.candidates=0, candidates="", delta.1.2=0)[0,]
    for (recruit in recruits) { 
        row <- Ljk.r[[label]][recruit,]
        which.L <- do.call(method,list(row))
        L <- min(row[which.L])  # a priori the low negLogL
        matches <- mothers[which.L]
        mom.L <- if (drop.Inf && all(matches == Inf))
                     "no match L==Inf" 
                 else {
                     if (return.full.match || length(matches) < .max.matches.to.report)
                         paste(collapse=",", matches)
                     else paste("***",length(matches),"matches")
                 }
        ans <- rbind(ans, data.frame(recruit=recruit, n.mother=length(matches), mother=mom.L, 
                                     data.source=data.source, negLogL=L, 
                                     candidate.window="", n.candidates=0, candidates="", delta.1.2=0))
    }
    rownames(ans) <- ans$recruit
    ans
}


extract.via.L.default <- function(Ljk.r, 
                                  data.source=.data.sources,
                                  return.full.match=FALSE,
                                  return.candidates=FALSE,
                                  drop.Inf=TRUE)
{

    # This default method allows us to produce reasonable diagnostics because
    # we know the way we are extracting candidates.  Returns a data frame with
    # recruits extracted from Ljk.r according to their negLogLik, as specified
    # by the method function.  If there is more than one mother, the mothers
    # are separated by comma.

    data.source <- match.arg(data.source)
    recruits <- dimnames(Ljk.r[[1]])[[1]]
    mothers <- dimnames(Ljk.r[[1]])[[2]]
    label <- switch(data.source, pericarp.priority = "Ljk.recruit.priority",
                                 use.all = "Ljk.recruit.use.all",
                                 pericarp = "Ljk.pericarp",
                                 seedling = "Ljk.seedling")
    ans <- data.frame(recruit="", n.mother=0, mother="", data.source="", negLogL=0, 
                      candidate.window="", n.candidates=0, candidates="", delta.1.2=0)[0,]
    for (recruit in recruits) { 
        ans <- rbind(ans, 
                     extract.via.L.default.recruit(recruit, 
                                                   Ljk.r, 
                                                   data.source=data.source,
                                                   return.full.match=return.full.match,
                                                   return.candidates=return.candidates,
                                                   drop.Inf=drop.Inf)
                    )
    }
    rownames(ans) <- ans$recruit
    ans
}


##################

extract.via.L.default.recruit <- function(recruit,
                                          Ljk.r, 
                                          data.source=.data.sources,
                                          return.full.match=FALSE,
                                          return.candidates=FALSE,
                                          drop.Inf=TRUE)
{
    # This default method allows us to produce reasonable diagnostics because
    # we know the way we are extracting candidates.  Returns a data frame with
    # recruits extracted from Ljk.r according to their negLogLik, as specified
    # by the method function.  If there is more than one mother, the mothers
    # are separated by comma.

    data.source <- match.arg(data.source)
    #recruits <- dimnames(Ljk.r[[1]])[[1]]
    mothers <- dimnames(Ljk.r[[1]])[[2]]
    label <- switch(data.source, pericarp.priority = "Ljk.recruit.priority",
                                 use.all = "Ljk.recruit.use.all",
                                 pericarp = "Ljk.pericarp",
                                 seedling = "Ljk.seedling")
    row <- Ljk.r[[label]][recruit,]
    L <- min(row)  # the low negLogL
    matches <- mothers[row == L]  # mothers matching the low negLogL
    mom.L <- if (drop.Inf && all(matches == Inf))  # the string representing the matching mothers
                    "no match L==Inf" 
             else {
                    if (return.full.match || length(matches) < .max.matches.to.report)
                        paste(collapse=",", matches)
                    else paste("***",length(matches),"matches")
             }

    candidates <- n.candidates <- c()

    # create backward probabilities, calculate LOD scores from those
    erow <- 10^row; min.erow <- min(erow)
    lod.scores <- log10(erow / min.erow)
    for (L.window in .included.candidate.window) {
        if (L == Inf) {
            n.candidates <- c(n.candidates, 0)
            candidates <- c(candidates, "")
        } else {
###HERE
        }
        which.win <- row >= L & row <= (L + L.window)
        n.candidates <- sum(which.win)
        candidates <- mothers[which.win]
    }
    # now calculate delta
    L.vals <- unique(sort(row))
    if (length(L.vals) == 1)
        delta.1.2 <- 0
    else delta.1.2 <- L.vals[1] - L.vals[2]
    ans <- list(recruit=recruit, n.mother=length(matches), mother=mom.L, 
                data.source=data.source, negLogL=L, 
                candidate.window=paste(collapse=",",round(.included.candidate.window,4)),
                n.candidates=paste(collapse=",",n.candidates), 
                candidates=if (return.candidates) paste(collapse=",",candidates) else "", 
                delta.1.2=0)
    ####
    ####
    ans
}


extract.via.L.old <- function(Ljk.r, data.source=.data.sources,
                                 method="min",
                                 return.full.match=FALSE,
                                 drop.Inf=TRUE)
{
    # Returns a data frame with recruit<tab>mother<tab>L.
    # If there is more than one mother, the mothers are separated by comma.

    data.source <- match.arg(data.source)
    recruits <- dimnames(Ljk.r[[1]])[[1]]
    mothers <- dimnames(Ljk.r[[1]])[[2]]
    label <- switch(data.source, pericarp.priority = "Ljk.recruit.priority",
                                 use.all = "Ljk.recruit.use.all",
                                 pericarp = "Ljk.pericarp",
                                 seedling = "Ljk.seedling")
    ans <- data.frame(recruit="",n.mother=0,mother="",data.source="",negLogL=0)[0,]
    for (recruit in recruits) { 
        row <- Ljk.r[[label]][recruit,]
        L <- do.call(method,list(row))
        matches <- mothers[which(row == L)]
        mom.L <- if (drop.Inf && (L == Inf || L == -Inf))
                     "no match L==Inf" 
                 else {
                     if (return.full.match || length(matches) < .max.matches.to.report)
                         paste(collapse=",", matches)
                     else paste("***",length(matches),"matches")
                 }
        ans <- rbind(ans, data.frame(recruit=recruit, n.mother=length(matches), mother=mom.L, 
                                     data.source=data.source, negLogL=L))
    }
    rownames(ans) <- ans$recruit
    ans
}


