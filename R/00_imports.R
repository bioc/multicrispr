#' @importFrom  assertive  assert_all_are_dirs
#' @importFrom  assertive  assert_all_are_existing_files
#' @importFrom  assertive  assert_all_are_false
#' @importFrom  assertive  assert_all_are_less_than
#' @importFrom  assertive  assert_all_are_matching_regex  assert_all_are_true
#' @importFrom  assertive  assert_are_identical  assert_are_same_length
#' @importFrom  assertive  has_names   assert_has_names
#' @importFrom  assertive  is_identical_to_true  is_scalar  is_windows
#' @importFrom  assertive  assert_is_subset  assert_has_no_duplicates
#' @importFrom  assertive  assert_is_all_of    assert_is_any_of  
#' @importFrom  assertive  assert_is_a_bool    assert_is_a_number
#' @importFrom  assertive  assert_is_a_string  assert_is_character
#' @importFrom  assertive  assert_is_numeric   assert_is_scalar
#' @importFrom  Biostrings            complement  DNAStringSet   
#' @importFrom  Biostrings            vmatchPattern  vcountPDict
#' @importFrom  BiocGenerics          end  end<-  invertStrand
#' @importFrom  BiocGenerics          start  start<-  strand  strand<-
#' @importFrom  BSgenome              getSeq   getBSgenome
#' @importFrom  data.table            fwrite
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
#' @importFrom  methods               as
#' @importFrom  tidyr                 separate_rows
#' @importFrom  utils                 download.file   getFromNamespace  head  
#' @importFrom  utils                 tail  read.csv  read.table
#' @importFrom  stats                 complete.cases
#' @importFrom  stringi               stri_detect_regex  stri_locate_all_fixed
#' @importFrom  stringi               stri_locate_all_regex   
#' @importFrom  stringi               stri_replace_first_fixed
#' @importFrom  stringi               stri_startswith_fixed
#' @importFrom  tidyselect            starts_with
NULL
