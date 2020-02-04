#' Extract subranges
#' 
#' Extract subranges from a \code{\link{GRanges-class}} object
#' 
#' Extract subranges (defined \code{\link{IRanges-class}} format) from a
#' \code{\link{GRanges}} object
#' 
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @param ir \code{\link[IRanges]{IRanges-class}}: subranges to be extracted
#' @param plot TRUE or FALSE (default)
#' @return \code{\link[GenomicRanges]{GRanges-class}}. 
#' @examples 
#' gr <- GenomicRanges::GRanges(c(A = 'chr1:1-100:+', B = 'chr1:1-100:-'))
#' gr$targetname <- 'AB'
#' ir <- IRanges::IRanges(c(A = '1-10', A = '11-20', B = '1-10', B = '11-20'))
#' extract_subranges(gr, ir, plot = TRUE)
#' @export
extract_subranges <- function(gr, ir, plot = FALSE){
  
    # Comply / Assert
    substart <- subwidth <- NULL
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(ir, 'IRanges')
    assert_has_names(gr)
    assert_has_names(ir)
    assert_is_subset(unique(names(ir)), names(gr))
    assert_is_a_bool(plot)
    
    # Convert and merge
    idt <- data.table::as.data.table(ir)
    gdt <- data.table::as.data.table(gr) %>% cbind(names = names(gr))
    gdt$width <- NULL
    setnames(idt, c(   'start',    'end',    'width'), 
                  c('substart', 'subend', 'subwidth'))
    mdt <- merge(gdt, idt, by = 'names')

    # Extract
    mdt[strand=='+', start := start-1+substart]
    mdt[strand=='+', end   := start-1+subwidth]
    mdt[strand=='-', end   := end+1-substart]
    mdt[strand=='-', start := end+1-subwidth]
    mdt[, c('substart', 'subend', 'subwidth') := NULL]
    mr <- dt2gr(mdt, seqinfo = seqinfo(gr))
    names(mr) %<>% uniquify()
    
    # Plot    
    if (plot)  print(plot_intervals(mr))
    
    # Return
    mr
}

#' Extract matching subranges
#' 
#' Extract subranges that match pattern
#' 
#' @param gr       \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome \code{\link{BSgenome}{BSgenome-class}}
#' @param pattern  string: search pattern in extended IUPAC alphabet
#' @param plot     TRUE or FALSE (default)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # PE example
#' #------------
#' require(magrittr)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                         HBB  = 'chr11:5227002:-',             # snp
#'                         HEXA = 'chr15:72346580-72346583:-',   # del
#'                         CFTR = 'chr7:117559593-117559595:+'), # ins
#'                       bsgenome)
#' gr %<>% extend_for_pe()
#' pattern <- strrep('N',20) %>% paste0('NGG')
#' extract_matchranges(gr, bsgenome, pattern, plot = TRUE)
#' 
#' # TFBS examples
#' #--------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#' gr <- bed_to_granges(bedfile, 'mm10') %>% extend()
#' extract_matchranges(gr, bsgenome, pattern = strrep('N',20) %>% paste0('NGG'))
#' @export
extract_matchranges <- function(gr, bsgenome, pattern, plot = FALSE){
  
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_has_names(gr)
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_a_string(pattern)

    # Extract
    matches <- unlist(vmatchPattern(pattern, getSeq(bsgenome, gr), fixed=FALSE))
    extract_subranges(gr, matches, plot = plot) %>% sort(ignore.strand = TRUE)
    
}


#' Find crispr spacers in targetranges
#' @param gr         \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param spacer     string: spacer pattern in extended IUPAC alphabet
#' @param pam        string: pam pattern in extended IUPAC alphabet
#' @param complement TRUE (default) or FALSE: also search in compl ranges?
#' @param plot       TRUE (default) or FALSE
#' @param verbose    TRUE (default) or FALSE
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#'     plot_intervals(gr)
#'     find_pe_spacers(gr, bsgenome)
#'     find_spacers(extend_for_pe(gr), bsgenome, complement = FALSE)
#'           # complement = FALSE because extend_for_pe  already 
#'           # adds  reverse complements and does so in a strand-specific 
#'           # manner
#'     
#' # TFBS example
#' #-------------
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10') %>% extend()
#'     find_spacers(gr, bsgenome)
#' @seealso \code{\link{find_pe_spacers}} to find prime editing spacers 
#' @export 
find_spacers <- function(
    gr, bsgenome, spacer = strrep('N', 20), pam = 'NGG', complement = TRUE, 
    verbose = TRUE, plot = TRUE
){

    if (complement){
        gr %<>% add_inverse_strand(plot = FALSE, verbose = verbose)
    }
    sites <- extract_matchranges(gr, bsgenome, paste0(spacer, pam)) %>% unique()
    spacers <-     extend(sites,  0, -3, bsgenome = bsgenome)
    pams    <- down_flank(sites, -2,  0, bsgenome = bsgenome)
    spacers$crisprname   <- names(spacers)
    spacers$crisprspacer <- BSgenome::getSeq(bsgenome, spacers, as.character=TRUE)
    spacers$crisprpam    <- BSgenome::getSeq(bsgenome, pams,    as.character=TRUE)
    spacers %>% sort(ignore.strand = TRUE)
    if (plot){
        print(plot_intervals(spacers, y='crisprname'))
        spacers$sitename <- NULL
    }
    spacers
}

#' Extend ranges for prime editing
#' 
#' Extend target ranges to span in which to look for spacer-pam seqs
#' 
#' Extend target ranges to find nearby spacers for prime editing
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param nrt       number: reverse transcription length
#' @param spacer    string: spacer pattern in extended IUPAC alphabet
#' @param pam       string: pam pattern in extended IUPAC alphabet
#' @param plot      TRUE (default) or FALSE
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome = bsgenome)
#'     find_pe_spacers(gr, bsgenome)
#'     (grext <- extend_for_pe(gr))
#'     find_spacers(grext, bsgenome, complement = FALSE)
#' @export
extend_for_pe <- function(
    gr, 
    bsgenome, 
    nrt    = 16, 
    spacer = strrep('N', 20), 
    pam    = 'NGG', 
    plot   = TRUE
){
    fw <- copy( gr, start=end(gr)+1-nrt-17, end=start(gr)-1+6,      strand='+')
    rv <- copy( gr, start=end(gr)+1-6,      end=start(gr)-1+nrt+17, strand='-')
    names(fw) %<>% paste0('_f')
    names(rv) %<>% paste0('_r')
    ext <- c(fw, rv)
    if (plot){
      gr$set <- 'PE target'
      fw$set <- "potential '+' spacers"
      rv$set <- "potential '-' spacers"
      print(plot_intervals(c(fw, gr, rv), color_var = 'set', y = 'set'))
    }
    ext
}


# Convert PAM into regex format
# @examples
# x <- 'NGG'
# pam_to_regex <- function(x){
#     assert_is_a_string(x)
#     x %>% stringi::stri_replace_all_fixed('N', '[ACGT]')
# }


# This regex-based function fails to find all crispr spacers
# It is, therefore, no longer used, just kept for reference purposes.
# Example explanation of why it fails
#     x <- Biostrings::DNAString('AGCAGCTGGGGCAGTGGTGGGGGGCCTTGGCGGCTACA')
#     stringi::stri_locate_all_regex(as.character(x), '[ACGT]{21}GG')
#     Biostrings::matchPattern('NNNNNNNNNNNNNNNNNNNNNGG', x, fixed = FALSE)
# regex_based_find_crispr_spacers <- function(
#   gr, bsgenome, pam = 'NGG', plot = TRUE, verbose = TRUE){
# 
#     # Assert. Import. Comply
#     assert_is_all_of(gr, 'GRanges')
#     assert_is_all_of(bsgenome, 'BSgenome')
#     assert_is_a_bool(verbose)
#     start <- substart <- crispr_start <- NULL
#     end <- subend <- crispr_end <- strand <- seqnames <- NULL
#     gr %<>% add_seq(bsgenome)
#     gr %<>% name_uniquely()
# 
#     # Find crispr sites in targetranges
#     pattern <- paste0('[ACGT]{20}', pam_to_regex(pam))
#     if (verbose) cmessage('\tFind %s crispr sites', pattern)
#     targetdt <- as.data.table(gr)
#     res <- targetdt$seq %>% stri_locate_all_regex(pattern)
#     cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
#     cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
#     targetdt [ , substart := vapply( res, cextract1, character(1)) ]
#     targetdt [ , subend   := vapply( res, cextract2, character(1)) ]
# 
#     # Rm crispr-free targetranges
#     idx <- targetdt[, substart == 'NA']
#     if (sum(idx)>0){
#         if (verbose)  cmessage('\t\tRm %d ranges with no crispr site', sum(idx))
#         targetdt %<>% extract(!idx)
#     }
# 
#     # Transform into crispr ranges
#     targetdt[, targetstart  := start]
#     targetdt[, targetend    := end  ]
# 
#     sites_dt  <-  tidyr::separate_rows(targetdt, substart, subend)          %>%
#                 data.table()                                                %>%
#                 extract(, substart := as.numeric(substart))                 %>%
#                 extract(, subend   := as.numeric(subend))                   %>%
#                 extract(, seq      := substr(seq, substart, subend))        %>%
#                 extract( strand=='+', crispr_start := start + substart - 1) %>%
#                 extract( strand=='+', crispr_end   := start + subend   - 1) %>%
#                 extract( strand=='-', crispr_start := end   - subend   + 1) %>%
#                 extract( strand=='-', crispr_end   := end   - substart + 1) %>%
#                 extract(, list(
#                             seqnames = seqnames, start = crispr_start,
#                             end = crispr_end,  strand = strand,  seq = seq,
#                             targetname = targetname, targetstart = targetstart,
#                             targetend = targetend))
#     sites <- GRanges(unique(sites_dt), seqinfo =  seqinfo(gr))
#     sites$site <- uniquify(sites$targetname)
#     spacer <- sites %>%     extend( 0, -3, stranded=TRUE, bsgenome = bsgenome)
#     pam    <- sites %>% down_flank(-2,  0, stranded=TRUE, bsgenome = bsgenome)
#     spacer$spacer <- spacer$seq
#     spacer$seq <- NULL
#     spacer$pam <- pam$seq
# 
#     # Plot. Message. Return
#     if (plot){
#         original  <- copy(sites, start=sites$targetstart, end=sites$targetend)
#         original$set <- 'target'
#         spacer$set <- 'spacer'
#         plot_intervals(c(original, spacer), color_var = 'set',
#                        size_var = 'set', y = 'site')
#         sites$set <- NULL
#     }
#     if (verbose)   cmessage('\t\t%d cas9 spacers across %d ranges',
#                         length(unique(spacer$spacer)), length(spacer))
#     return(sites)
# }



# Add three prime extensions
# 
# Add three prime extensions for prime editing
# 
# A prime editor contains spacer and 3'extensions
# This function adds these extension seqs for the specified spacer(range)s
# 
#   ====== spacer =====>.........................
#   1--------------17----GG--------->           .
#   <--------------------CC----------           .
#       <======  extension ==========............
#       <====pbs=====<===rttemplate==
#                    
# @param spacers   \code{\link[GenomicRanges]{GRanges-class}}: spacers
# @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
# @param fixes     character vector: '+' strand fix seqs
# @param nprimer   n primer nucleotides (default 13, max 17)
# @param nrt       n rev transcr nucleotides (default 16, recomm. 10-16)
# @param spacer    string: spacer pattern in extended IUPAC alphabet
# @param pam       string: pam pattern in extended IUPAC alphabet
# @param plot      TRUE (default) or FALSE
# @return  Prime editing spacer \code{\link[GenomicRanges]{GRanges-class}}, 
# with 3' extension sequences (required for cloning) added.
# @examples
# require(magrittr)
# bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
# gr <- GenomicRanges::GRanges(
#         seqnames = c(PRNP = 'chr20:4699600',             # snp
#                      HBB  = 'chr11:5227002',             # snp
#                      HEXA = 'chr15:72346580-72346583',   # del
#                      CFTR = 'chr7:117559593-117559595'), # ins
#         strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#         seqinfo  = BSgenome::seqinfo(bsgenome))
#         
# spacers <- gr %>%  extend_for_pe() %>% find_spacers(bsgenome)
#         
# BSgenome::getSeq(bsgenome, gr)
# add_prime_extensions(gr, bsgenome)
# @seealso \code{\link{find_spacers}} to find standard crispr sites
# @export
# add_prime_extensions <- function(
#   spacers, 
#   bsgenome, 
#   rttemplate = getSeq(bsgenome, down_flank(spacers, -2, -3+nrt), as.character = TRUE), 
#   nprimer = 13, 
#   nrt = 16, 
#   spacer =  strrep('N', 20), 
#   pam = 'NGG', 
#   plot = TRUE
# ){
#     
#     # Assert
#     assert_is_all_of(spacers, 'GRanges')
#     assert_is_all_of(bsgenome, 'BSgenome')
#     assert_is_character(fixes)
#     assert_all_are_matching_regex(fixes, '^[ACGTacgt]+$')
#     assert_is_a_number(nprimer)
#     assert_is_a_number(nrt)
#     assert_all_are_less_than(nprimer, 17)
# .

#     # Extensions
#     extension <- invertStrand(down_flank(spacers, -2-nprimer, -3+nrt))
#     extension$seq  <- get_plus_seq(bsgenome, extension)              # Get "+" seq
#     substr( extension$seq, 
#             extension$tstart-start(extension)+1, 
#             extension$tend  -start(extension)+1) <- fixes[extension$tname]
#     extension$seq[as.logical(strand(extension)=='-')] %<>% revcomp() # Revcomp for "-" seqs
# 
#     # Plot
#     if (plot){
#         spacer$part<-'spacer'; primer$part<-'primer'; rtranscripts$part<-'rtranscripts'
#         allranges <- c(spacer, primer, rtranscripts)
#         allranges$part %<>% factor(rev(c('spacer', 'rtranscripts', 'primer')))
#         plot_intervals(allranges, y = 'site', color_var = 'part', 
#             size_var = 'part', facet_var = c('seqnames', 'targetname'))
#         spacer$part <- primer$part <- rtranscripts$part <- NULL
#     }
#     
#     # Add sequences and return
#     names(mcols(spacer)) %<>% stri_replace_first_fixed('seq', 'spacer')
#     spacer$pam <- pam$seq 
#     spacer$extension <- extension$seq
#     spacer
# }


