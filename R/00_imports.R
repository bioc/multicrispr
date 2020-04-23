#' @importFrom  assertive.base        assert_all_are_true  is_identical_to_true
#' @importFrom  assertive.base        assert_all_are_false
#' @importFrom  assertive.files       assert_all_are_existing_files
#' @importFrom  assertive.numbers     assert_all_are_less_than
#' @importFrom  assertive.properties  has_names   assert_has_names
#' @importFrom  assertive.reflection  is_windows
#' @importFrom  assertive.sets        assert_is_subset
#' @importFrom  assertive.strings     assert_all_are_matching_regex
#' @importFrom  assertive.types       assert_is_all_of    assert_is_any_of  
#' @importFrom  assertive.types       assert_is_a_bool    assert_is_a_number
#' @importFrom  assertive.types       assert_is_a_string  assert_is_character
#' @importFrom  Biostrings            complement  DNAStringSet   
#' @importFrom  Biostrings            vmatchPattern  vcountPDict
#' @importFrom  BiocGenerics          end  end<-  invertStrand
#' @importFrom  BiocGenerics          start  start<-  strand  strand<-
#' @importFrom  BSgenome              getSeq   getBSgenome
#' @importFrom  data.table            :=  data.table  as.data.table  setnames
#' @importFrom  data.table            setnames  setorderv   setnafill    .SD
#' @importFrom  GenomeInfoDb          genome
#' @importFrom  GenomeInfoDb          seqinfo   seqinfo<-
#' @importFrom  GenomeInfoDb          seqlevels  seqlevels<-  seqlevelsInUse
#' @importFrom  GenomeInfoDb          seqlevelsStyle  seqlevelsStyle<- 
#' @importFrom  GenomeInfoDb          seqnames  seqnames<-   
#' @importFrom  GenomeInfoDb          sortSeqlevels  standardChromosomes
#' @importFrom  GenomicRanges         GRanges  mcols  mcols<-
#' @importFrom  ggplot2               aes   aes_string
#' @importFrom  ggplot2               facet_wrap  geom_point  geom_segment  
#' @importFrom  ggplot2               ggplot  ggtitle  scale_x_continuous
#' @importFrom  ggplot2               scale_alpha_manual  scale_size_manual
#' @importFrom  ggplot2               theme_bw  unit  xlab  ylab
#' @importFrom  grid                  arrow
#' @importFrom  magrittr              %>%   %<>%   and
#' @importFrom  magrittr              extract  extract2  set_names
#' @importFrom  methods               as
#' @importFrom  tidyr                 separate_rows
#' @importFrom  utils                 download.file   getFromNamespace  head  
#' @importFrom  utils                 tail  read.csv  read.table
#' @importFrom  stringi               stri_detect_regex  stri_locate_all_fixed
#' @importFrom  stringi               stri_locate_all_regex   
#' @importFrom  stringi               stri_replace_first_fixed
#' @importFrom  stringi               stri_startswith_fixed
NULL
