% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IC.R
\name{IC}
\alias{IC}
\title{IC object}
\usage{
## Create IC object
ICobj <- makeIC(RAobj, samID=NULL,
                filter=list(MAF=0.05, MISS=0.2, BIN=100, DEPTH=5, PVALUE=0.01, MAXDEPTH=500))

## Functions (Methods) of an IC object
ICobj$computeMap(chrom = NULL, init_r = 0.01, ep = 0.001, method = "EM", err = T, mapped = T, nThreads = 1)
ICobj$createLG(LODthres = 10, nComp = 10)
ICobj$maskSNP(snps)
ICobj$mergeLG(LG)
ICobj$orderLG(chrom = NULL, mapfun = "haldane", weight = "LOD2", ndim = 30, spar = NULL)
ICobj$plotChr(mat = "rf", filename = NULL, chrS = 2, lmai = 2)
ICobj$plotLG(LG = NULL, mat = "rf", filename = NULL, interactive = TRUE)
ICobj$plotLM(LG = NULL, fun = "haldane", col = "black")
ICobj$plotSyn()
ICobj$print(what = NULL, ...)
ICobj$removeLG(LG)
ICobj$removeSNP(snps)
ICobj$rf_2pt(nClust = 2, err=FALSE)
ICobj$unmaskSNP(snps)
ICobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)
}
\description{
Class for storing RA data and associated functions for analysis of intercross (F2) populations.
}
\details{
An IC object is created from the \code{\link{makeIC}} function and contains RA data,
various statistics of the dataset that have been computed, and functions (methods)
for analyzing the data. Information in an IC object are specific to intercross (F2) populations.
}
\section{Methods(Functions)}{

\describe{
\item{\code{\link{$createLG}}}{Create linkage group(s).}
\item{\code{\link{$computeMap}}}{Compute linkage maps for a given marker order.}
\item{\code{\link{$maskSNP}}}{Mask SNP(s) in the dataset.}
\item{\code{\link{$mergeLG}}}{Merge linkage groups.}
\item{\code{\link{$orderLG}}}{Order the SNPs in the linkage group(s).}
\item{\code{\link{$plotChr}}}{Plot the heatmap of the 2-point recombination fraction 
estimates (or LOD scores) when SNPs are ordered according to the genome assembly.}
\item{\code{\link{$plotLG}}}{Plot the heatmap of the 2-point recombination fraction 
estimates (or LOD scores) when SNPs are ordered according to their linkage groups.}
\item{\code{\link{$plotLM}}}{Plot linkage maps.}
\item{\code{\link{$plotSyn}}}{Produce a Synteny plot.}
\item{\code{\link{$print}}}{Print summary information for FS object.}
\item{\code{\link{$removeLG}}}{Remove linkage group(s).}
\item{\code{\link{$removeSNP}}}{Remove SNP(s) from linkage group(s).}
\item{\code{\link{$rf_2pt}}}{Compute the 2-point recombination fraction (and LOD score) between all SNP pairs.}
\item{\code{\link{$unmaskSNP}}}{Unmask SNP(s) in the data set.}
\item{\code{\link{$writeLM}}}{Write linkage mapping results to a file.}
}
}

\author{
Timothy P. Bilton
}
\section{Super class}{
\code{\link[GUSbase:RA]{GUSbase::RA}} -> \code{IC}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-IC-new}{\code{IC$new()}}
\item \href{#method-IC-print}{\code{IC$print()}}
\item \href{#method-IC-removeSNP}{\code{IC$removeSNP()}}
\item \href{#method-IC-removeLG}{\code{IC$removeLG()}}
\item \href{#method-IC-mergeLG}{\code{IC$mergeLG()}}
\item \href{#method-IC-maskSNP}{\code{IC$maskSNP()}}
\item \href{#method-IC-unmaskSNP}{\code{IC$unmaskSNP()}}
\item \href{#method-IC-rf_2pt}{\code{IC$rf_2pt()}}
\item \href{#method-IC-createLG}{\code{IC$createLG()}}
\item \href{#method-IC-orderLG}{\code{IC$orderLG()}}
\item \href{#method-IC-plotLG}{\code{IC$plotLG()}}
\item \href{#method-IC-plotChr}{\code{IC$plotChr()}}
\item \href{#method-IC-plotLM}{\code{IC$plotLM()}}
\item \href{#method-IC-plotSyn}{\code{IC$plotSyn()}}
\item \href{#method-IC-computeMap}{\code{IC$computeMap()}}
\item \href{#method-IC-writeLM}{\code{IC$writeLM()}}
\item \href{#method-IC-cometPlot}{\code{IC$cometPlot()}}
\item \href{#method-IC-rocketPlot}{\code{IC$rocketPlot()}}
\item \href{#method-IC-RDDPlot}{\code{IC$RDDPlot()}}
\item \href{#method-IC-clone}{\code{IC$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="extractVar"><a href='../../GUSbase/html/RA.html#method-RA-extractVar'><code>GUSbase::RA$extractVar()</code></a></span></li>
<li><span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="mergeSamples"><a href='../../GUSbase/html/RA.html#method-RA-mergeSamples'><code>GUSbase::RA$mergeSamples()</code></a></span></li>
<li><span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="writeRA"><a href='../../GUSbase/html/RA.html#method-RA-writeRA'><code>GUSbase::RA$writeRA()</code></a></span></li>
<li><span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="writeVCF"><a href='../../GUSbase/html/RA.html#method-RA-writeVCF'><code>GUSbase::RA$writeVCF()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-new"></a>}}
\if{latex}{\out{\hypertarget{method-IC-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$new(R6obj)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-print"></a>}}
\if{latex}{\out{\hypertarget{method-IC-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$print(what = NULL, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-removeSNP"></a>}}
\if{latex}{\out{\hypertarget{method-IC-removeSNP}{}}}
\subsection{Method \code{removeSNP()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$removeSNP(snps)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-removeLG"></a>}}
\if{latex}{\out{\hypertarget{method-IC-removeLG}{}}}
\subsection{Method \code{removeLG()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$removeLG(LG)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-mergeLG"></a>}}
\if{latex}{\out{\hypertarget{method-IC-mergeLG}{}}}
\subsection{Method \code{mergeLG()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$mergeLG(LG)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-maskSNP"></a>}}
\if{latex}{\out{\hypertarget{method-IC-maskSNP}{}}}
\subsection{Method \code{maskSNP()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$maskSNP(snps)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-unmaskSNP"></a>}}
\if{latex}{\out{\hypertarget{method-IC-unmaskSNP}{}}}
\subsection{Method \code{unmaskSNP()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$unmaskSNP(snps)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-rf_2pt"></a>}}
\if{latex}{\out{\hypertarget{method-IC-rf_2pt}{}}}
\subsection{Method \code{rf_2pt()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$rf_2pt(nClust = 2, err = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-createLG"></a>}}
\if{latex}{\out{\hypertarget{method-IC-createLG}{}}}
\subsection{Method \code{createLG()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$createLG(LODthres = 10, nComp = 10, reset = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-orderLG"></a>}}
\if{latex}{\out{\hypertarget{method-IC-orderLG}{}}}
\subsection{Method \code{orderLG()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$orderLG(
  LG = NULL,
  mapfun = "haldane",
  weight = "LOD2",
  ndim = 30,
  spar = NULL,
  filename = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-plotLG"></a>}}
\if{latex}{\out{\hypertarget{method-IC-plotLG}{}}}
\subsection{Method \code{plotLG()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$plotLG(LG = NULL, mat = "rf", filename = NULL, interactive = TRUE, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-plotChr"></a>}}
\if{latex}{\out{\hypertarget{method-IC-plotChr}{}}}
\subsection{Method \code{plotChr()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$plotChr(chrom = NULL, mat = c("rf"), filename = NULL, chrS = 2, lmai = 2)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-plotLM"></a>}}
\if{latex}{\out{\hypertarget{method-IC-plotLM}{}}}
\subsection{Method \code{plotLM()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$plotLM(LG = NULL, fun = "haldane", col = "black")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-plotSyn"></a>}}
\if{latex}{\out{\hypertarget{method-IC-plotSyn}{}}}
\subsection{Method \code{plotSyn()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$plotSyn()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-computeMap"></a>}}
\if{latex}{\out{\hypertarget{method-IC-computeMap}{}}}
\subsection{Method \code{computeMap()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$computeMap(
  chrom = NULL,
  init_r = 0.001,
  ep = 0.001,
  method = "EM",
  err = TRUE,
  multiErr = FALSE,
  mapped = TRUE,
  nThreads = 1,
  inferOPGP = TRUE,
  rfthres = 0.1,
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-writeLM"></a>}}
\if{latex}{\out{\hypertarget{method-IC-writeLM}{}}}
\subsection{Method \code{writeLM()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-cometPlot"></a>}}
\if{latex}{\out{\hypertarget{method-IC-cometPlot}{}}}
\subsection{Method \code{cometPlot()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$cometPlot(
  filename = NULL,
  cex = 1,
  maxdepth = 500,
  maxSNPs = 1e+05,
  res = 300,
  color = NULL,
  ind = FALSE,
  ncores = 1,
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-rocketPlot"></a>}}
\if{latex}{\out{\hypertarget{method-IC-rocketPlot}{}}}
\subsection{Method \code{rocketPlot()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$rocketPlot(
  ploid = 2,
  filename = NULL,
  cex = 1,
  maxdepth = 500,
  maxSNPs = 1e+05,
  res = 300,
  scaled = TRUE,
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-RDDPlot"></a>}}
\if{latex}{\out{\hypertarget{method-IC-RDDPlot}{}}}
\subsection{Method \code{RDDPlot()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$RDDPlot(filename = NULL, maxdepth = 500, maxSNPs = 1e+05, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-IC-clone"></a>}}
\if{latex}{\out{\hypertarget{method-IC-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IC$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
