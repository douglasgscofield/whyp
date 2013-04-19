
#### #### #### #### #### #### #### ####
Pr.recruit.at.random <- function(recruit, dat, method=c("smouse","new"), Pr.null=.Pr.null)
{
    method <- match.arg(method)
    loci <- attr(dat$iseedling, "locus.names")
    n.loci <- attr(dat$iseedling, "n.loci")
    locus.columns <- attr(dat$iseedling, "locus.columns")
    if (method == "smouse") {
        Pr.p <- Pr.pericarp.at.random(recruit, dat, Pr.null=Pr.null)
        Pr.s <- Pr.seedling.at.random(recruit, dat, method="smouse")
        if (length(Pr.p) != length(Pr.s)) 
            stop("length of pericarp and seedling vectors must be the same")
        # where we have a pericarp genotype (not 1), make the seedling not
        # matter (=1)
        Pr.s[Pr.p < 1] <- 1
        Pr.random <- Pr.s * Pr.p
    } else if (method == "new") {
        missing.val <- -1
        Pr.p <- Pr.pericarp.at.random(recruit, dat, Pr.null=Pr.null, 
                                      missing.val=missing.val)
        Pr.s <- Pr.seedling.at.random(recruit, dat, method="new")
        if (length(Pr.p) != length(Pr.s)) 
            stop("length of pericarp and seedling vectors must be the same")
        # where we have a pericarp genotype (not 1), make the seedling not
        # matter (=1)
        Pr.s[Pr.p != missing.val] <- 1
        Pr.p[Pr.p == missing.val] <- 1
        Pr.random <- Pr.s * Pr.p
    }
    ####
    Pr.random
}

Pr.seedling.at.random <- function(recruit, dat, method=c("smouse","new"),
                                  missing.val=1)
{
    method <- match.arg(method)
    loci <- attr(dat$iseedling, "locus.names")
    n.loci <- attr(dat$iseedling, "n.loci")
    locus.columns <- attr(dat$iseedling, "locus.columns")
    Pr.random <- rep(0, n.loci)
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]; l2 <- l1 + 1
        s1 <- dat$iseedling[recruit,l1]; s2 <- dat$iseedling[recruit,l2]
        if (s1 == missing.allele.index || s2 == missing.allele.index) 
            ans <- missing.val
        else {
            if (method == "smouse") {
                if (s1 == s2) {  # homozygote: freq of s1 in population
                    ans <- dat$allele.table[loci[locus],s1] *
                           dat$allele.table[loci[locus],s1]
                } else {  # heterozygote: 1/2 freq of each allele in population
                    ans <- 2.0 * dat$allele.table[loci[locus],s1] * 
                           dat$allele.table[loci[locus],s2]
                }
            } else if (method == "new") {
                if (s1 == s2) {  # homozygote: freq of s1 in population
                    ans <- dat$allele.table[loci[locus],s1]
                } else {  # heterozygote: 1/2 freq of each allele in population
                    ans <- 0.5 * (dat$allele.table[loci[locus],s1] + 
                            dat$allele.table[loci[locus],s2])
                }
            }
        }
        Pr.random[locus] <- ans
    }
    ####
    Pr.random
}

Pr.pericarp.at.random <- function(recruit, dat, Pr.null=.Pr.null, missing.val=1)
{
    loci <- attr(dat$ipericarp, "locus.names")
    n.loci <- attr(dat$ipericarp, "n.loci")
    locus.columns <- attr(dat$iseedling, "locus.columns")
    Pr.random <- rep(0, n.loci)
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]; l2 <- l1 + 1
        s1 <- dat$ipericarp[recruit,l1]; s2 <- dat$ipericarp[recruit,l2]
        if (s1 == missing.allele.index || s2 == missing.allele.index) 
            ans <- missing.val
        else {
            if (s1 == s2) { # homozygote: (1-Pr.null)*pi*pi + Pr.null*pi*(1-pi)
                ans <- ((1 - Pr.null) * dat$allele.table[loci[locus],s1] * 
                      dat$allele.table[loci[locus],s1]) +
                     (Pr.null * dat$allele.table[loci[locus],s1] *
                     (1 - dat$allele.table[loci[locus],s1]))
            } else {
                ans <- 2* dat$allele.table[loci[locus],s1] * 
                        dat$allele.table[loci[locus],s2]
            }
        }
        Pr.random[locus] <- ans
    }
    ####
    Pr.random
}
#### #### #### #### #### #### #### ####

# functions for calculating backward probabilities


Pr.back.mother.given.recruit <- function(mother, recruit, dat)
# produce backward probability values for mother | seedling
{
    this.iseedling <- dat$iseedling[recruit,]
    this.ipericarp <- dat$ipericarp[recruit,]
    this.imother <- dat$imother[mother,]
    allele.table <- dat$allele.table
    loci <- attr(dat$iseedling, "locus.names"); n.loci <- attr(dat$iseedling, "n.loci")
    locus.columns <- attr(dat$iseedling, "locus.columns")
    seedling.loci.missing <- pericarp.loci.missing <- rep(FALSE, n.loci)
    names(seedling.loci.missing) <- names(pericarp.loci.missing) <- loci
    Pr.back <- rep(0, n.loci); names(Pr.back) <- loci
    for (locus in 1:n.loci) {
        l1 <- locus.columns[locus]; l2 <- l1 + 1
        m1 <- this.imother[l1]; m2 <- this.imother[l2]
        p1 <- this.ipericarp[l1]; p2 <- this.ipericarp[l2]
        s1 <- this.iseedling[l1]; s2 <- this.iseedling[l2]
        Pr.back[locus] <- Pr.back.locus.mother.given.recruit(m1, m2, s1, s2, p1, p2, loci[locus], allele.table)
    }
    ####
    ####
    prod(Pr.back)
}


Pr.back.locus.mother.given.recruit <- function(m1, m2, s1, s2, p1, p2, locus.name, allele.table)
{
    Pr.null <- .Pr.null
    mother.missing <- if (missing.allele.index %in% c(m1,m2)) TRUE else FALSE
    Pr.par.seed <- Pr.back.locus.parent.given.seedling(m1, m2, s1, s2, locus.name, allele.table)
    if (mother.missing) {
        Prob <- 1  # could be anything
    } else if (missing.allele.index == p1 || missing.allele.index == p2) {
        # missing data in pericarp, only based on seedling
        Prob <- Pr.par.seed
    } else {
        Prob <- 0
        if (s1 == s2) { # seedling homozygote
            if (p1 == p2) { # pericarp homozygote
                if (m1 == m2) { # mother homozygote
                    if (m1 != s1) { # homozygotes don't match
                        Prob <- 0  
                    } else { # homozygotes match
                        Prob <- 1 - Pr.null*(1 - Pr.par.seed)
                    }
                } else { # mother heterozygote
                    if (m1 == s1 || m2 == s1) { # one mother allele matches,
                        # could be null genotyping error
                        Prob <- Pr.null*(1 - Pr.par.seed)
                    } else { # no chance mother made seedling
                        Prob <- 0  
                    }
                }
            } else { # pericarp heterozygote
                if ((p1 == m1 && p2 == m2) || (p1 == m2 && p2 == m1)) {
                    # het pericarp and mother match, no genotyping error assumed
                    Prob <- 1
                } else {
                    Prob <- 0
                }
            }
        } else { # seedling heterozygote
            if (p1 == p2) { # pericarp homozygote
                # reorder seedling alleles so first is the one that matches pericarp
                if (p1 == s2) { tmp <- s2; s2 <- s1; s1 <- tmp }
                # backward probabilities
                Pr.hh.hi <- Pr.back.locus.parent.given.seedling(s1, s1, s1, s2, locus.name, allele.table)
                Pr.hi.hi <- Pr.back.locus.parent.given.seedling(s1, s2, s1, s2, locus.name, allele.table)
                Pr.hx.hi <- (1 - Pr.hh.hi - Pr.hi.hi)
                if (m1 == m2) { # mother homozygote
                    if (m1 == p1) { # homozygotes match
                        Prob <- 1 - Pr.null*(1 - Pr.hh.hi)
                    } else { # homozygotes do not match, could be homozygote
                        # of other seedling allele, or some other allele, either
                        # way it can't match with our error scheme
                        Prob <- 0
                    }
                } else { # mother heterozygote
                    if ((p1 == m1 && p2 == m2) || (p1 == m2 && p2 == m1)) { # mother
                        # heterozygote matches seedling heterozygote
                        Prob <- Pr.null*(1 - Pr.hh.hi)*(Pr.hi.hi/(Pr.hi.hi + Pr.hx.hi))
                    } else if (s1 == m1) { # mother heterozygote and shares the shared allele
                        # between seedling and pericarp, but not the other
                        Prob <- Pr.null*(1 - Pr.hh.hi)*(Pr.hx.hi/(Pr.hi.hi + Pr.hx.hi))
                    } else { # mother heterozygote but caries the unshared allele
                        Prob <- 0
                    }
                }
            } else { # pericarp heterozygote
                if ((p1 == m1 && p2 == m2) || (p1 == m2 && p2 == m1)) {
                    # het pericarp and mother match, no genotyping error assumed
                    Prob <- 1
                } else {
                    Prob <- 0
                }
            }
        }
    }
    ####
    ####
    Prob
}


Pr.back.locus.parent.given.seedling <- function(par1, par2, seed1, seed2, locus.name, allele.table)
# produce backward probability per-locus values for mother | seedling
{
    Prob <- 0
    parent.missing <- if (missing.allele.index %in% c(par1,par2)) TRUE else FALSE
    if (parent.missing) {
        Prob <- 1  # could be anything
    } else if (missing.allele.index == seed1 || missing.allele.index == seed2) {
        # Hardy-Weinberg proportion based on parent alleles
        p <- (allele.table[locus.name,par1] * allele.table[locus.name,par2])
        Prob <- if (par1 != par2) 2*p else p
    } else {
        Prob <- 0
        if (seed1 == seed2) { # seedling homozygote
            if (par1 == par2) { # mother homozygote
                if (par1 != seed1) { # homozygotes don't match
                    Prob <- 0  
                } else { # homozygotes match
                    Prob <- (allele.table[locus.name,seed1] *
                                allele.table[locus.name,seed1] *
                                allele.table[locus.name,seed1])
                }
            } else { # mother heterozygote
                if (par1 == seed1 || par2 == seed1) { # one mother allele matches,
                    # father provided the other
                    Prob <- (allele.table[locus.name,seed1] * 
                                allele.table[locus.name,seed1] * 
                                (1 - allele.table[locus.name,seed1]))
                } else { # no chance mother made seedling
                    Prob <- 0  
                }
            }
        } else { # seedling heterozygote
            if (par1 == par2) { # mother homozygote
                if (par1 == seed1) { # father provided seed2
                    Prob <- (allele.table[locus.name,par1] *
                                allele.table[locus.name,par1] *
                                allele.table[locus.name,seed2])
                } else if (par1 == seed2) { # father provided seed1
                    Prob <- (allele.table[locus.name,par1] *
                                allele.table[locus.name,par1] *
                                allele.table[locus.name,seed1])
                } else { # no chance mother made seedling
                    Prob <- 0  
                }
            } else { # mother heterozygote
                if ((par1 == seed1 && par2 == seed2) ||
                    (par1 == seed2 && par2 == seed1)) { # matching heterozygotes,
                    # father provided the other allele, don't know which
                    Prob <- (allele.table[locus.name,seed1] * 
                                allele.table[locus.name,seed2] *
                                (allele.table[locus.name,seed1] + 
                                allele.table[locus.name,seed2]))
                } else if (par1 == seed1 || par2 == seed1 ||
                            par1 == seed2 || par2 == seed2) { # father provided seed1 or seed2
                    Prob <- (allele.table[locus.name,seed1] * 
                                allele.table[locus.name,seed2] *
                                (1 - allele.table[locus.name,seed1] - 
                                    allele.table[locus.name,seed2]))
                } else { # no chance mother made seedling
                    Prob <- 0  
                }
            }
        }
    }
    ####
    ####
    Prob
}


