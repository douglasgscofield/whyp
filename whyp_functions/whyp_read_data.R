# functions for reading and preparing mixed assay data

read.mixedassay.data <- function(file, file.format=c("genalex"), data.source)
{
    cat("read.mixedassay.data:",data.source,"data, file",file,"format",file.format,"\n")
    file.format <- match.arg(file.format)
    data.source <- match.arg(data.source, c("mother","pericarp","seedling"))
    data <- switch(file.format,
        genalex = readGenalex(file)
    )
    names(data)[names(data) == attr(data,"sample.title")] <- "id"
    names(data)[names(data) == attr(data,"pop.title")] <- "data.source"
    data[["data.source"]] <- data.source
    attr(data,"pop.title") <- NULL
    attr(data,"data.source") <- data.source
    data
}


read.mother <- function(file=file.mother)
{ read.mixedassay.data(file=file, data.source="mother") }


read.pericarp <- function(file=file.pericarp)
{ read.mixedassay.data(file=file, data.source="pericarp") }


read.seedling <- function(file=file.seedling)
{ read.mixedassay.data(file=file, data.source="seedling") }


merge.recruit <- function(pericarp, seedling, 
                          mod.loci=TRUE,
                          missing.action=c("exclude","include"))
{
    # They could have different sets of loci, so re-format each to match a
    # merged format.  There could also be individuals in one dataset and not in
    # the other, so add entries containing missing data to those
    missing.action <- match.arg(missing.action)
    p.loci <- attr(pericarp,"locus.names")
    s.loci <- attr(seedling,"locus.names")
    r.loci <- sort(unique(c(p.loci,s.loci)))
    n.loci <- length(r.loci)
    ploidy <- attr(pericarp,"ploidy")
    p.new <- pericarp[,1:2]
    s.new <- seedling[,1:2]
    names(p.new)[1] <- names(s.new)[1] <- "id"
    if (length(unique(p.new$id)) < length(p.new$id))
        stop("duplicate IDs in pericarp data")
    if (length(unique(s.new$id)) < length(s.new$id))
        stop("duplicate IDs in seedling data")
    # add columns for loci that aren't there
    for (r.locus in r.loci) {
        if (r.locus %in% p.loci)
            p.new <- cbind(p.new, getGenalexLocus(pericarp, r.locus))
        else {  # empty locus
            loc <- list()
            loc[[r.locus]] <- loc[[paste(sep=".",r.locus,"2")]] <- 
                rep(missing.allele[1],nrow(p.loci))
            full.p <- cbind(full.p, as.data.frame(loc))
        }
        if (r.locus %in% s.loci)
            s.new <- cbind(s.new, getGenalexLocus(seedling, r.locus))
        else {  # empty locus
            loc <- list()
            loc[[r.locus]] <- loc[[paste(sep=".",r.locus,"2")]] <- 
                rep(missing.allele[1],nrow(r.loci))
            new.r <- cbind(new.r, as.data.frame(loc))
        }
    }
    ids <- sort(union(p.new$id, s.new$id))
    missing.p.ids <- ids[! ids %in% p.new$id]
    missing.s.ids <- ids[! ids %in% s.new$id]
    if (missing.action == "include") {
        # add rows for genotypes that aren't there
        if (length(missing.p.ids) > 0) {
            tmp <- matrix(missing.allele[1],nrow=length(missing.p.ids),ncol=n.loci*ploidy)
            tmp <- data.frame(id=missing.p.ids, data.source="pericarp", tmp)
            names(tmp) <- names(p.new)
            p.new <- rbind(p.new,tmp)
            p.new <- p.new[order(p.new$id),]
        }
        if (length(missing.s.ids) > 0) {
            tmp <- matrix(missing.allele[1],nrow=length(missing.s.ids),ncol=n.loci*ploidy)
            tmp <- data.frame(id=missing.s.ids, data.source="seedling", tmp)
            names(tmp) <- names(s.new)
            s.new <- rbind(s.new,tmp)
            s.new <- s.new[order(s.new$id),]
        }
        cat(missing.action,"non-matching entries in datasets: added",
            length(missing.p.ids),"blank entry(s) to pericarp,",
            length(missing.s.ids),"to seedling\ntotal",
            length(ids),"entries merged per dataset\n")
    } else if (missing.action == "exclude") {
        # remove IDs missing in one or the other data set
        ids <- sort(intersect(p.new$id, s.new$id))
        #ids <- sort(ids[! ids %in% unique(c(missing.p.ids, missing.s.ids))])
        if (length(ids) == 0) stop("no ids in common between the datasets")
        s.new <- subset(s.new, id %in% ids)[order(ids),]
        p.new <- subset(p.new, id %in% ids)[order(ids),]
        cat(missing.action,"non-matching entries in datasets: removed",
            length(missing.p.ids),"unique entry(s) from pericarp,",
            length(missing.s.ids),"from seedling\ntotal",
            length(ids),"entries kept per dataset\n")
    } else stop("unknown missing.action")
    attr(p.new,"n.loci") <- n.loci
    attr(p.new,"ploidy") <- ploidy
    attr(p.new,"locus.names") <- r.loci
    attr(p.new,"locus.columns") <- sort(which(names(p.new) %in% r.loci))
    attr(p.new,"n.samples") <- nrow(p.new)
    attr(p.new,"genetic.data.format") <- "genalex"
    rownames(p.new) <- 1:nrow(p.new)
    attributes(s.new) <- attributes(p.new)
    attr(p.new,"data.source") <- "pericarp"
    attr(s.new,"data.source") <- "seedling"
    list(pericarp=p.new, seedling=s.new)
}


index.genotypes <- function(data, allele.table, return.matrix=TRUE)
{
    # Among several things, make sure the columns match the order of the loci
    # in allele.table, and don't barf if the table includes a locus that isn't
    # in the data.
    ai <- attr(allele.table,"index.list")
    ids <- data[,1]
    tloci <- dimnames(allele.table)[[1]] # attr(data,"locus.names")
    loci <- attr(data, "locus.names")
    if (! all(loci %in% tloci))
        stop("Not all loci in dataset in the loci to index")
    # if loci in data are a proper subset of loci in the allele table, order
    # them as in the table, but don't include the ones that aren't in the table
    # (of course)
    data <- reorderGenalexLoci(data, loci)
    newdata <- data
    for (locus in loci) {
        ldata <- getGenalexLocus(data, locus)
        newldata <- matrix(missing.allele.index, nrow=nrow(ldata),ncol=ncol(ldata))
        for (row in 1:nrow(ldata)) {
            any.missing <- any(ldata[row,] %in% missing.allele)
            for (col in 1:ncol(ldata)) {
                allele <- ldata[row,col]
                newldata[row,col] <- 
                    if (any.missing && missing.single.allele.method == "as.missing" &&
                        ! allele %in% missing.allele) {
                        cat("recruit", data$id[row], "locus", locus, 
                            "non-missing allele made missing due to",
                            missing.single.allele.method,"method\n")
                        newldata[row,col] <- missing.allele.index
                    } else if (allele %in% missing.allele) {
                        newldata[row,col] <- missing.allele.index
                    } else { 
                        newldata[row,col] <- ai[[locus]][[as.character(allele)]] 
                    }
        }}
        newdata <- putGenalexLocus(newdata, locus, newldata)
    }
    if (return.matrix) {
        newdata <- as.matrix(newdata[,3:ncol(newdata)])
    }
    rownames(newdata) <- ids
    attr(newdata,"locus.names") <- loci
    attr(newdata,"n.loci") <- length(loci)
    attr(newdata,"locus.columns") <- locus.columns <- which(dimnames(newdata)[[2]] %in% loci)
    n.missing.loci <- rep(0,length(ids))
    name.missing.loci <- rep("",length(ids))
    for (r in seq(along=ids)) {
        tmp <- (newdata[r,locus.columns] == missing.allele.index)
        name.missing.loci[r] <- paste(sep=",",collapse=",",loci[tmp])
        n.missing.loci[r] <- sum(tmp)
    }
    attr(newdata,"n.missing.loci") <- n.missing.loci
    attr(newdata,"name.missing.loci") <- name.missing.loci
    newdata
}

