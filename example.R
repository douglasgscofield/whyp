options(prompt="whyp-examples> ", 
        stringsAsFactors=FALSE, 
        windowsBuffered=FALSE)
reload <- function(doit=FALSE) if (doit) source("examples.R")


source("whyp.R")
source("readGenalex.R")
source("whyp/whyp_alleles.R")
source("whyp/whyp_backward_probability.R")
source("whyp/whyp_check.R")
source("whyp/whyp_extract_likelihood.R")
source("whyp/whyp_likelihood.R")
source("whyp/whyp_read_data.R")
source("whyp/whyp_utility.R")


whyp.analysis <- function(dat) {

  ans <- whyp(dat)

  # Alternatives:
  #
  # recruits <- dat$recruit$pericarp$id[1:100]
  # ans <- whyp(dat, recruits)  # only do whyp analysis for the first 100 recruits
  
  # ans <- whyp(recruits, dat, heur.tissue="pericarp")  # assume mismatches may occur in the pericarp

  # ans <- subset(ans, method == "pericarp")  # only return the analyses based on pericarp tissue

  return(ans)
}

# Collect genotype data from adults and pericarps, produce reports on foreign alleles,
# genotype mismatches, and potential null alleles, and drop a problematic locus while
# loading data

assemble.v1 <- function()
  assemble.whyp.data(file.mother="genalex-format-adult-genotypes.txt",
                    file.pericarp="genalex-format-pericarp-genotypes.txt",
                    report.foreign.recruit="foreign-allele-report.txt",
                    report.mismatch.recruit="potential-mismatch-report.txt",
                    report.null.recruit="potential-null-allele-report.txt",
                    read.freqs=FALSE,
                    drop.loci=c("ZAG110")
                    )

# Collect genotype data from adults, pericarps and seedlings, produce reports on foreign alleles,
# genotype mismatches, and potential null alleles, and keep all loci.

assemble.v2 <- function()
  assemble.whyp.data(file.mother="genalex-format-adult-genotypes.txt",
                    file.pericarp="genalex-format-pericarp-genotypes.txt",
                    file.seedling="genalex-format-seedling-genotypes.txt",
                    report.foreign.recruit="foreign-allele-report.txt",
                    report.mismatch.recruit="potential-mismatch-report.txt",
                    report.null.recruit="potential-null-allele-report.txt",
                    read.freqs=FALSE,
                    drop.loci=NULL
                    )


