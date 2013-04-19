# utility routines


# print genotypes of recruit(s) and associated mother(s)
mr <- function (rr,m=NULL,dat,table=FALSE) 
{
  ans <- data.frame()
  for (r in rr) {
    ans <- rbind(ans, dat$recruit$pericarp[dat$recruit$pericarp[,1] == r,])
    ans <- rbind(ans, dat$recruit$seedling[dat$recruit$seedling[,1] == r,])
  }
  if (nrow(ans) == 0)
    stop("none of the recruits were found in the dataset")
  if (is.null(m)) # extract mom from recruit name(s)
    m <- unique(matrix(unlist(strsplit(rr,"_")), nrow=2)[1, ])  
  ans <- rbind(dat$mother[dat$mother[,1] %in% m, ], ans)
  if (table)
    write.table(ans, file="", quote=FALSE, row.names=FALSE, sep="\t")
  else 
    print(ans)
}



list.to.data.frame <- function(data, create.freqs=TRUE, cols=c("locus","allele","freq"))
{
    # creates a three-column data.frame
    if (length(cols) != 3) stop("cols must have 3 names")
    ans <- data.frame()
    by.1 <- names(data)
    for (b in by.1) {
        a <- unlist(data[[b]]); a <- a[order(names(a))]
        if (create.freqs) a <- a/sum(a)
        locdf <- list()
        locdf[[cols[1]]] <- b
        locdf[[cols[2]]] <- names(a)
        locdf[[cols[3]]] <- as.vector(a)
        ans <- rbind(ans, data.frame(locdf))
    }
    ans
}


data.frame.to.list <- function(data, by)
{
    # creates an order-two list.
    if (ncol(data) != 3) stop("data.frame must have 3 cols")
    if (length(by) != 2) stop("by must have 2 names")
    by.1 <- sort(unique(as.character(data[[by[1]]])))
    by.2.name <- split(data[[by[2]]], data[[by[1]]])
    by.2.val <- split(data[,3], data[[by[1]]])
    val.list <- list()
    for (item in by.1) {
        vals <- as.vector(by.2.val[[item]])
        names(vals) <- as.vector(by.2.name[[item]])
        val.list[[item]] <- as.list(vals)
        #missing.list[[item]] <- 0
    }
    val.list
}


