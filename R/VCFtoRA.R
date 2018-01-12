
VCFtoRA <- function(infilename){
  
  outfilename <- paste0(infilename,".ra.tab")
  
  headerlist = c('CHROM', 'POS')
  
  Lines <- readLines(infilename)
  
  for(i in 1:length(Lines)){
    line = Lines[[i]]
    line = trimws(line)
    
    
    if(substr(line,1,2) == "##")
      next
    else if(substr(line,1,6) == "#CHROM"){
      stop("test")
      ## extract the sample IDS
      line = strsplit(line, split="\t")
      headerlist = c(headerlist, line[[1]][10:length(line[[1]])])
      ## Write the header to the RA file
      write(paste(paste(headerlist, collapse="\t"),sep="\n"), file=outfilename, ncolumns = 1)
      cat("found",length(headerlist),"samples\n")
    }
      
    
  }
  
}

infilename = "//isamba/dataset/MBIE_genomics4production/active/manuka/tassel_ref/06_VCF/manuka_filtered_aln.vcf"
