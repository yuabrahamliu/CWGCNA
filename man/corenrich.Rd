% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emethyl_enrich.R
\name{corenrich}
\alias{corenrich}
\title{Perform correlation-based gene functional enrichment or EnrichR enrichment}
\usage{
corenrich(
  hyperprobes = NULL,
  hypoprobes = NULL,
  pairedRNA = NULL,
  pairedmethyl = NULL,
  abscut = 0.7,
  generegions = c("TSS200"),
  platform = 450,
  dbs = c("GO_Biological_Process_2018", "BioPlanet_2019", "Reactome_2016"),
  write = FALSE
)
}
\arguments{
\item{hyperprobes}{The vector with the hypermethylated methylation probes or
DOWN-regulated RNA genes to be analyzed. Default is NULL.}

\item{hypoprobes}{The vector with the hypomethylated methylation probes or
UP-regulated RNA genes to be analyzed. Default value is NULL. Can transfer
2 vectors to the 2 parameters \code{hyperprobes} and \code{hypoprobes}, or
transfer 1 vector to 1 of these 2 parameters and leave the other as NULL.}

\item{pairedRNA}{The RNA part of the paired RNA-DNA methylation dataset. If
both this parameter and \code{pairedmethyl} are not NULL, a correlation-
based method will be used to find the genes whose RNA expression level in
the RNA data significantly correlated with any of the \code{hyperprobes}
and \code{hypoprobes} in the paired methylation data. The probe values in
the methylation data will be converted to gene values first, and then the
correlation analysis will be performed. In addition to conducting EnrichR
enrichment for the correlated genes finally identified, each of the genes
will also get its own function annotated. Should transfer a matrix to this
parameter, with columns representing samples and rows representing genes.
Row names are gene symbols and column names are sample IDs. Default value
is NULL, and in this case, the correlation-based method will not be used,
and instead, function annotation and enrichment analysis will be directly
conducted on the hyper and hypomethylated methylation probes transferred
to \code{hyperprobes} and \code{hypoprobes}, with these probes mapped to
corresponding genes first. Or, if \code{hyperprobes} and \code{hypoprobes}
have already been gene symbols, the mapping step will be skipped.}

\item{pairedmethyl}{The methylation part of the paired RNA-DNA methylation
dataset. Similar to \code{pairedRNA}, if both these 2 parameters are not
NULL, the correlation-based method will be used, but if any one of them is
NULL, EnrichR enrichment will be directly conducted on probes or genes in
\code{hyperprobes} and \code{hypoprobes}. \code{pairedmethyl} should be a
beta value matrix with columns representing samples and rows for PROBEs.
Row names are PROBE names. Column names are sample names. The sample names
should be the same as the ones in \code{pairedRNA}, because they are data
for paired samples. Default value of \code{pairedmethyl} is NULL.}

\item{abscut}{When the correlation method is used to select the genes with a
significant correlation to the probes, the ones with an absolute Pearson
correlation coefficient value greater than this parameter value will be
selected and used for function analysis. Default is 0.7.}

\item{generegions}{When \code{hyperprobes} and \code{hypoprobes} contain the
names of methylation probes, this parameter will be used when converting
them into genes. That is, if a probe located in the regions defined by it,
the corresponding gene will be used for the function analysis. Default is
"TSS200", so probes within TSS200 regions will be used for gene mapping.
It can also be a vector, such as \code{c("TSS200", "TSS1500", "1stExon")},
so that probes within any of these 3 regions will be used for the mapping.
This gene mapping process will be used when any of \code{pairedRNA} and
\code{pairedmethyl} is NULL, so the EnrichR enrichment will be performed
directly on the genes mapped from the parameters \code{hyperprobes} and
\code{hypoprobes}. When both \code{pairedRNA} and \code{pairedmethyl} are
not NULL, the correlation-based analysis will be used, and this mapping
will also be used because before screening for the correlated genes, the
the methylation probes will be merged into gene methylation values first,
and then used for calculating the correlation with the gene RNA values in
the paired RNA data. The gene methylation values are summarized according
to the regions in this parameter \code{generegions}.}

\item{platform}{The platform of the probes. Can be 27 (for 27k platform),
450 (for 450k platform), or 850 (for EPIC platform). Default is 450.}

\item{dbs}{The databases for gene functional enrichment analysis. Default is
\code{c("GO_Biological_Process_2018", "BioPlanet_2019", "Reactome_2016")}.}

\item{write}{A logical value indicating whether the gene function results
need to be written into txt files in the working directory. The default
value is FALSE.}
}
\value{
A list with four slots recording the gene annotation and functional
enrichment results for the enhanced and inhibited genes (when correlation-
based method is used) or for the genes mapped from the hypomethylated and
hypermethylated probes (when correlation-based method is not used). When
\code{hyperprobes} and \code{hypoprobes} contain gene names rather than
probes, the slot \code{hypergeneenrich} records the result of the genes in
\code{hyperprobes}, which are DOWN-regulated genes. In constrast, the slot
\code{hypogeneenrich} recordes the results for \code{hypoprobes}, which
UP-regulated genes. If \code{write} is TRUE, txt files will be generated
to save these results.
}
\description{
Perform correlation-based gene functional enrichment or EnrichR enrichment.
}
