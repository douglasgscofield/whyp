# functions for checking for nulls and mismatches

check.recruit.null <- function(dat=NULL, mother=NULL, pericarp=NULL, report=FALSE)
{
    pericarp <- if (! is.null(dat)) dat$recruit$pericarp 
                else if (! is.null(pericarp)) pericarp 
                     else stop("pericarp missing")
    mother <- if (! is.null(dat)) dat$mother
                else if (! is.null(mother)) mother 
                     else stop("mother missing")
    if (attr(pericarp,"ploidy") != 2) stop("ploidy must be 2")
    loci <- attr(pericarp,"locus.names"); n.loci <- attr(pericarp,"n.loci")
    locus.columns <- attr(pericarp, "locus.columns")
    pericarp.names <- pericarp[,1]; n.pericarp <- length(pericarp.names)
    mother.names <- mother[,1]; n.mother <- length(mother.names)

    ans <- list()
    ans$n.potential.null.loci <- ans$name.potential.null.loci <- array(0, dim=c(n.pericarp,n.mother), dimnames=list(pericarp.names,mother.names))
    ans$name.potential.null.loci[] <- ""

    for (j in 1:n.pericarp) {
        pericarp.printed <- FALSE
        for (k in 1:n.mother) {
            potential.null.loci <- c()
            match.pmatch <- 0
            for (locus in 1:n.loci) {
                l1 = locus.columns[locus]; l2 = l1 + 1
                m1 = mother[k,l1]; m2 = mother[k,l2]
                p1 = pericarp[j,l1]; p2 = pericarp[j,l2]

                # See if we match at this locus.  If we do, continue.  If we
                # don't, then check for potential null in the pericarp, which I
                # assume can occur when the pericarp is a homozygote, the mom
                # is a heterozygote, and the pericarp matches one of the mom
                # alleles.  If that's not true, we definitely don't match, so
                # stop comparing to this mother.

                if (all(sort(c(m1,m2)) == sort(c(p1,p2)))) {
                    match.pmatch <- match.pmatch + 1
                    next
                } else if (p1 == p2 && m1 != m2 && (p1 == m1 || p1 == m2)) {
                    match.pmatch <- match.pmatch + 1
                    potential.null.loci <- c(potential.null.loci, loci[locus])
                } else break
            }
            if (match.pmatch == n.loci) {
                ans$n.potential.null.loci[j,k] <- length(potential.null.loci)
                ans$name.potential.null.loci[j,k] <- paste(collapse=",", 
                                                           potential.null.loci)
                if (report && length(potential.null.loci) > 0) {
                    if (! pericarp.printed) {
                        pericarp.printed <- TRUE
                        printGenalexGenotype(pericarp, j, 
                                               callout.locus=potential.null.loci, 
                                               sep="\t")
                    }
                    printGenalexGenotype(mother, k, 
                                           callout.locus=potential.null.loci, 
                                           sep="\t")
                }
            }
        }
        if (report && sum(ans$n.potential.null.loci[j,])) cat("\n")
        if ((j %% 50) == 0) cat("finished null check through pericarp",j,"...\n")
    }
    cat("finished null check for", n.pericarp, "pericarps\n")
    for (a in c("n.potential.null.loci", "name.potential.null.loci")) {
        attr(dat$recruit$pericarp,a) <- attr(dat$ipericarp,a) <- ans[[a]]
    }
    #attr(dat$recruit$pericarp,"n.potential.null.loci") <- attr(dat$ipericarp,"n.potential.null.loci") <- n.potential.null.loci
    #attr(dat$recruit$pericarp,"name.potential.null.loci") <- attr(dat$ipericarp,"name.potential.null.loci") <- name.potential.null.loci
    ####
    ####
    invisible(dat)
}


report.recruit.null <- function(dat, pericarp.name=NULL, pericarp.idx=NULL, 
                                report.full=FALSE, sep="\t")
{
    p.idx <- if (is.null(pericarp.idx)) 
                 which(dat$recruit$pericarp[,1] %in% pericarp.name)
             else pericarp.idx
    n.mat <- attr(dat$recruit$pericarp,"n.potential.null.loci")
    name.mat <- attr(dat$recruit$pericarp,"name.potential.null.loci")
    if (report.full) 
        cat(paste(collapse=sep, names(dat$recruit$pericarp)[c(1:2,attr(dat$recruit$pericarp,"locus.columns"))]),"\n")
    for (p in p.idx) {
        null.moms <- which(n.mat[p,] > 0)
        if (length(null.moms)) {
            printGenalexGenotype(dat$recruit$pericarp, p, sep=sep)
            for (mom in null.moms) {
                null.names <- strsplit(name.mat[p, mom],",")[[1]]
                printGenalexGenotype(dat$mother, mom, callout.locus=null.names, 
                                    sep=sep)
            }
            if (report.full) cat("\n")
        }
    }
}


edit.recruit.mismatch <- function(dat, method=c("set.all.missing",
                                                "set.pericarp.missing",
                                                "set.seedling.missing"),
                                  update.attributes=TRUE)
{
    method <- match.arg(method)
    if (is.null(attr(dat$ipericarp,"n.mismatch.loci")))
        stop("n.mismatch.loci not present, must run check.recruit.mismatch first")
    n.mm <- attr(dat$ipericarp,"n.mismatch.loci")
    name.mm <- attr(dat$ipericarp,"name.mismatch.loci")
    p.loci <- attr(dat$ipericarp,"locus.columns")
    names(p.loci) <- attr(dat$ipericarp,"locus.names")
    s.loci <- attr(dat$iseedling,"locus.columns")
    names(s.loci) <- attr(dat$iseedling,"locus.names")
    ploidy <- attr(dat$recruit$pericarp,"ploidy")
    m.ats <- list(); m.ats$p <- m.ats$s <- list()
    for (a in c("n.missing.loci", "name.missing.loci")) {
        m.ats$p[[a]] <- attr(dat$ipericarp,a)
        m.ats$s[[a]] <- attr(dat$iseedling,a)
    }
    for (a in c("n.double.missing.loci", "name.double.missing.loci")) {
        m.ats[[a]] <- attr(dat$ipericarp,a)
    }

    for (r in 1:nrow(dat$ipericarp)) {
        if (n.mm[r] == 0) next
        if (name.mm[r] == "") stop("n.mismatch.loci > 0 but no locus names")

        mm.loci <- unlist(strsplit(name.mm[r],","))

        if (method == "set.pericarp.missing" || method == "set.all.missing") {
            cols <- p.loci[mm.loci]
            if (ploidy > 1) for (a in 1:(ploidy - 1)) cols <- c((cols + a), cols)
            #if (any(dat$ipericarp[r,cols] %in% missing.allele.index))
            #    stop("missing data in mismatched ipericarp locus, row ",r)
            dat$ipericarp[r,cols] <- missing.allele.index
            # update missing-data attributes
            m.ats$p$n.missing.loci[r] <- m.ats$p$n.missing.loci[r] + length(mm.loci)
            m.ats$p$name.missing.loci[r] <- if (m.ats$p$name.missing.loci[r] != "") 
                paste(collapse=",", c(m.ats$p$name.missing.loci[r], mm.loci)) 
            else paste(collapse=",", mm.loci)
        }
        if (method == "set.seedling.missing" || method == "set.all.missing") {
            cols <- s.loci[mm.loci]
            if (ploidy > 1) for (a in 1:(ploidy - 1)) cols <- c(cols, (cols + a))
            #if (any(dat$iseedling[r,cols] %in% missing.allele.index))
            #    stop("missing data in mismatched iseedling locus, row ",r)
            dat$iseedling[r,cols] <- missing.allele.index
            # update missing-data attributes
            m.ats$s$n.missing.loci[r] <- m.ats$s$n.missing.loci[r] + length(mm.loci)
            m.ats$s$name.missing.loci[r] <- if (m.ats$s$name.missing.loci[r] != "") 
                paste(collapse=",", c(m.ats$s$name.missing.loci[r], mm.loci)) 
            else paste(collapse=",", mm.loci)
        }
        if (method == "set.all.missing") {
            # update *.double.missing.loci attributes here
            m.ats$n.double.missing.loci[r] <- m.ats$n.double.missing.loci[r] + length(mm.loci)
            m.ats$name.double.missing.loci[r] <- if (m.ats$name.double.missing.loci[r] != "") 
                paste(collapse=",", c(m.ats$name.double.missing.loci[r], mm.loci)) 
            else paste(collapse=",", mm.loci)
        }
    }

    if (update.attributes) {
        for (a in c("n.missing.loci", "name.missing.loci")) {
            attr(dat$recruit$pericarp,a) <- attr(dat$ipericarp,a) <- m.ats$p[[a]]
            attr(dat$recruit$seedling,a) <- attr(dat$iseedling,a) <- m.ats$s[[a]]
        }
        for (a in c("n.double.missing.loci", "name.double.missing.loci")) {
            attr(dat$recruit$pericarp,a) <- attr(dat$recruit$seedling,a) <- attr(dat$ipericarp,a) <- attr(dat$iseedling,a) <- m.ats[[a]]
        }
    }
    attr(dat,"edit.recruit.mismatch") <- method
    attr(dat,"edit.recruit.mismatch.update.attributes") <- update.attributes
    dat
}


check.recruit.mismatch <- function(dat, report=FALSE, sep="\t")
{
    if (! all(dat$recruit$pericarp[,1] == dat$recruit$seedling[,1]))
        stop("pericarp & seedling entries out of sync, merge data first")
    # form intersection of locus names, compare only those
    p.loci <- attr(dat$recruit$pericarp,"locus.columns")
    names(p.loci) <- attr(dat$recruit$pericarp,"locus.names")
    s.loci <- attr(dat$recruit$seedling,"locus.columns")
    names(s.loci) <- attr(dat$recruit$seedling,"locus.names")
    locus.names <- intersect(names(p.loci), names(s.loci))
    ploidy <- attr(dat$recruit$pericarp,"ploidy")
    if (ploidy != 2) stop("can't handle ploidy != 2")
    p.loci <- p.loci[locus.names]; s.loci <- s.loci[locus.names]

    ans <- list()
    ans$n.mismatch.loci <- ans$n.double.missing.loci <- rep(0,nrow(dat$ipericarp))
    ans$name.mismatch.loci <- ans$name.double.missing.loci <- rep("",nrow(dat$ipericarp))

    if (report)
        cat(paste(collapse=sep, names(dat$recruit$pericarp)[c(1:2,attr(dat$recruit$pericarp,"locus.columns"))]),"\n")

    for (r in 1:nrow(dat$ipericarp)) {

        mismatch.loci <- double.missing.loci <- c()
        for (locus in locus.names) {

            # genotype data for this locus
            p.locus <- dat$recruit$pericarp[r,p.loci[locus]:(p.loci[locus]+ploidy-1)]
            s.locus <- dat$recruit$seedling[r,s.loci[locus]:(s.loci[locus]+ploidy-1)]

            # Check for a mismatch between pericarp and seedling genotypes.
            # This occurs when the seedling does not contain a pericarp allele,
            # and is not (according to the data) related to the pericarp at
            # this locus.
            if (! any(p.locus %in% missing.allele) &&
                ! any(s.locus %in% missing.allele) &&
                ! any(s.locus %in% p.locus))
                mismatch.loci <- c(mismatch.loci, locus)

            # if missing data in both pericarp and seedling
            if (all(p.locus %in% missing.allele) && 
                all(s.locus %in% missing.allele))
                double.missing.loci <- c(double.missing.loci, locus)

        }

        if (length(mismatch.loci) > 0) {
            ans$n.mismatch.loci[r] <- length(mismatch.loci)
            ans$name.mismatch.loci[r] <- paste(collapse=",",mismatch.loci)
            ans$n.double.missing.loci[r] <- length(double.missing.loci)
            ans$name.double.missing.loci[r] <- paste(collapse=",",double.missing.loci)
            if (report) {
                printGenalexGenotype(dat$recruit$pericarp, r, 
                                       callout.locus=mismatch.loci, sep=sep)
                printGenalexGenotype(dat$recruit$seedling, r, 
                                       callout.locus=mismatch.loci, sep=sep)
                cat("\n")
            }
        }
    }

    for (a in c("n.mismatch.loci", "name.mismatch.loci", 
                "n.double.missing.loci", "name.double.missing.loci")) {
        attr(dat$recruit$pericarp,a) <- attr(dat$recruit$seedling,a) <- attr(dat$ipericarp,a) <- attr(dat$iseedling,a) <- ans[[a]]
    }
    ####
    invisible(dat)
}

