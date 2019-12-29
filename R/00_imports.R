#' @importFrom  assertive.base        assert_all_are_true  is_identical_to_true
#' @importFrom  assertive.base        assert_all_are_false
#' @importFrom  assertive.files       assert_all_are_existing_files
#' @importFrom  assertive.numbers     assert_all_are_less_than
#' @importFrom  assertive.properties  has_names
#' @importFrom  assertive.sets        assert_is_subset
#' @importFrom  assertive.strings     assert_all_are_matching_regex
#' @importFrom  assertive.types       assert_is_all_of    assert_is_any_of  
#' @importFrom  assertive.types       assert_is_a_bool    assert_is_a_number
#' @importFrom  assertive.types       assert_is_a_string  assert_is_character
#' @importFrom  Biostrings            complement  DNAStringSet   vmatchPattern
#' @importFrom  BiocGenerics          end  end<-  invertStrand
#' @importFrom  BiocGenerics          start  start<-  strand  strand<-
#' @importFrom  BSgenome              getSeq
#' @importFrom  data.table            :=  data.table  as.data.table  setnames
#' @importFrom  data.table            setnames  setorderv
#' @importFrom  GenomeInfoDb          genome
#' @importFrom  GenomeInfoDb          seqinfo   seqinfo<-
#' @importFrom  GenomeInfoDb          seqlevels  seqlevels<-  seqlevelsInUse
#' @importFrom  GenomeInfoDb          seqlevelsStyle  seqlevelsStyle<- 
#' @importFrom  GenomeInfoDb          seqnames  seqnames<-   
#' @importFrom  GenomeInfoDb          sortSeqlevels  standardChromosomes
#' @importFrom  GenomicRanges         GRanges  mcols  mcols<-
#' @importFrom  ggplot2               aes_string  facet_wrap  geom_segment  
#' @importFrom  ggplot2               ggplot  ggtitle  theme_bw  unit  xlab  
#' @importFrom  ggplot2               ylab
#' @importFrom  grid                  arrow
#' @importFrom  magrittr              %>%   %<>%   extract  extract2  set_names
#' @importFrom  methods               as
#' @importFrom  tidyr                 separate_rows
#' @importFrom  utils                 getFromNamespace  head  tail
#' @importFrom  utils                 read.csv  read.table
#' @importFrom  stringi               stri_detect_regex  stri_locate_all_fixed
#' @importFrom  stringi               stri_locate_all_regex   
#' @importFrom  stringi               stri_replace_first_fixed
NULL
