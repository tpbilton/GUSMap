#' @export RA

### R6 class for data aligned to reference assembly
RA <- R6Class("RA",
              public = list(
                initialize = function(List){
                  private$genon     <- List$genon
                  private$ref       <- List$ref
                  private$alt       <- List$alt
                  private$chrom     <- List$chrom
                  private$pos       <- List$pos
                  private$SNP_Names <- List$SNP_Names
                  private$indID     <- List$indID
                  private$nSnps     <- List$nSnps
                  private$nInd      <- List$nInd
                  private$gform     <- List$gform
                  private$AFrq      <- List$AFrq
                },
                #### Diagonostic functions ####
                ## Ratio of alleles for heterozygous genotype calls (observed vs expected)
                cometPlot = function(model="random", alpha=NULL, filename=NULL, cex=1, maxdepth=500, ...){
                  cometPlot(private$ref, private$alt, model=model, alpha=alpha, filename=filename, cex=cex, maxdepth=maxdepth, ...)
                }
                ###############################
              ),
              private = list(
                genon = NULL,
                ref = NULL,
                alt = NULL,
                chrom = NULL,
                pos = NULL,
                SNP_Names = NULL,
                indID = NULL,
                nSnps = NULL,
                nInd = NULL,
                gform = NULL,
                AFrq = NULL
              )
)
