#' @importFrom  Biostrings            complement  DNAStringSet   
#' @importFrom  Biostrings            vmatchPattern  vcountPDict
#' @importFrom  Biostrings            writeXStringSet
#' @importFrom  BiocGenerics          end  end<-  invertStrand
#' @importFrom  BiocGenerics          start  start<-  strand  strand<-
#' @importFrom  BSgenome              getSeq   getBSgenome
#' @importFrom  data.table            fread   fwrite
#' @importFrom  data.table            :=  data.table  as.data.table  setnames
#' @importFrom  data.table            setnames  setorderv   setnafill    .SD
#' @importFrom  GenomeInfoDb          genome
#' @importFrom  GenomeInfoDb          seqinfo   seqinfo<-
#' @importFrom  GenomeInfoDb          seqlevels  seqlevels<-  seqlevelsInUse
#' @importFrom  GenomeInfoDb          seqlevelsStyle  seqlevelsStyle<- 
#' @importFrom  GenomeInfoDb          seqnames  seqnames<-   
#' @importFrom  GenomeInfoDb          sortSeqlevels  standardChromosomes
#' @importFrom  GenomicRanges         granges  GRanges  mcols  mcols<-
#' @importFrom  ggplot2               aes   aes_string
#' @importFrom  ggplot2               facet_wrap  geom_point  geom_segment  
#' @importFrom  ggplot2               ggplot  ggtitle  scale_x_continuous
#' @importFrom  ggplot2               scale_alpha_manual  scale_size_manual
#' @importFrom  ggplot2               scale_linetype_manual  scale_color_manual
#' @importFrom  ggplot2               theme_bw  unit  xlab  ylab
#' @importFrom  grid                  arrow
#' @importFrom  magrittr              %>%   %<>%   add   and
#' @importFrom  magrittr              extract  extract2  set_names
#' @importFrom  methods               as  is
#' @importFrom  tidyr                 separate_rows
#' @importFrom  utils                 download.file   getFromNamespace  head  
#' @importFrom  utils                 tail  read.csv  read.table
#' @importFrom  reticulate            py_module_available
#' @importFrom  stats                 complete.cases
#' @importFrom  stringi               stri_detect_fixed      stri_detect_regex
#' @importFrom  stringi               stri_locate_all_fixed  stri_locate_all_regex   
#' @importFrom  stringi               stri_replace_first_fixed
#' @importFrom  stringi               stri_startswith_fixed
#' @importFrom  tidyselect            starts_with
NULL
