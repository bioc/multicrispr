report_counts <- function(hits){
    
    if (length(hits)==0){
        cmessage('\t\t%s No cas9seqs which hit any target')
        return(NULL)
    }
    
    thits <- table(hits)
    
    most  <- thits[[ 1             ]]
    least <- thits[[ length(thits) ]]
    moststring  <-  most  %>% as.character()
    leaststring <-  least %>% 
                    stringi::stri_pad_left(' ', width = ceiling(log10(most)))
    
    # right tail, since distribution is geometric-like
    tailstring  <-  names(thits)[ length(thits) ]
    modestring  <-  names(thits)[1] %>% 
                    stringi::stri_pad_left(
                        ' ', 
                        width = ceiling(log10(as.numeric(tailstring))))  
                            # https://stackoverflow.com/a/47190862
    tail <- tailstring %>% as.numeric()
    mode <- modestring %>% as.numeric()
    
    cmessage('\t\t%s cas9seq%s hit%s %s target%s', 
            moststring,  
            ifelse(most ==1, ' ', 's'), 
            ifelse(most ==1, 's', ' ') , 
            modestring, 
            ifelse(mode==1, ' ', 's'))
    cmessage('\t\t%s cas9seq%s hit%s %s target%s', 
            leaststring, 
            ifelse(least==1, ' ', 's'), 
            ifelse(least==1, 's', ' ') , 
            tailstring, 
            ifelse(tail==1, ' ', 's'))
}

#' Count ontargets
#' @param cas9ranges   data.table(chr, start, end, strand)
#' @param targetseqs   character
#' @param mismatch     numeric(1)
#' @param bsgenome     BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param verbose      logical(1)
#' @return data.table
#' @examples 
#' require(magrittr)
#' bsgenome <-  BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' 
#' targetdt <- bedfile %>% read_bed() %>% flank_fourways(bsgenome)
#' cas9dt   <- bedfile %>% find_flankonly_cas9s(bsgenome)
#' cas9seqs   <- cas9dt$cas9seq %>% head(10)
#' 
#' cas9dt[cas9seqs, on = 'cas9seq'][, .N, by = 'cas9seq'][['N']]
#' count_ontargets(cas9seqs, targetseqs, bsgenome, mismatch=0)
#' cas9dt[cas9seqs[5], on='cas9seq'][, .N, by = 'cas9seq']
#' count_ontargets(cas9seqs[5], targetseqs, bsgenome, mismatch=0)
#' 
#' count_ontargets(cas9seqs, targetseqs, bsgenome, mismatch=1)
#' count_ontargets(cas9seqs, targetseqs, bsgenome, mismatch=2)
#' @export
count_ontargets <- function(
    cas9seqs, 
    targetseqs, 
    bsgenome,
    mismatch, 
    verbose = TRUE
){
    # Assert
    assertive.types::assert_is_character(cas9seqs)
    assertive.types::assert_is_character(targetseqs)
    assertive.types::assert_is_a_number(mismatch)
    assertive.types::assert_is_a_bool(verbose)
    
    # targetseqs
    cas9seqdt <- data.table::data.table(cas9seq = cas9seqs) %>% 
                 extract(, .(ncas9= .N), by = 'cas9seq')
    targetseqdt <-  data.table::data.table(targetseq = targetseqs) %>% 
                    extract(, .(ntarget = .N), by = 'targetseq')
    if (verbose)  cmessage(
                    '\tScan cas9seq ontargets with %d mismatches', mismatch)
    vres <- Biostrings::vcountPDict(
                Biostrings::DNAStringSet(cas9seqdt$cas9seq),
                Biostrings::DNAStringSet(targetseqdt$targetseq),
                min.mismatch = mismatch,
                max.mismatch = mismatch)
    
    vres %<>% set_rownames(cas9seqdt$cas9seq)
    vres %<>% set_colnames(targetseqdt$targetseq)
    dres <- vres %>%  
            data.table::data.table(keep.rownames = 'cas9seq') %>%
            data.table::melt(id.vars = 'cas9seq', variable.name = 'targetseq', value.name = 'npertarget') %>%
            merge(targetseqdt, by = 'targetseq') %>% 
            extract( , .(ontargets = sum(npertarget * ntarget)), by = 'cas9seq')
    
    if (verbose) report_counts(dres$ontargets)
    
    # cater for possibly duplicated cas9seqs
    hits <- dres[cas9seqs, ontargets, on = 'cas9seq']
    
    return(hits)
}





#' Scan genome
#' @param cas9seqs character(1)
#' @param bsgenome BSgenome object
#' @return data.table
#' @examples
#' require(magrittr)
#' bsgenome <-  BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' cas9seqs <- read_bed(bedfile)        %>%
#'             head(30)                 %>%  
#'             slop_fourways(bsgenome)  %>% 
#'             find_cas9s(bsgenome)     %>% 
#'             extract2('cas9seq')
#' @export
scan_genome <- function(
    cas9seqs, 
    bsgenome,
    mismatch, 
    verbose = TRUE
){
    # Assert
    assertive.base::assert_is_identical_to_true(
        methods::is(bsgenome, "BSgenome"))

    # Scan
    cas9seqdt <- data.table::data.table(
                    cas9seq = cas9seqs, 
                    index   = 1:length(cas9seqs))

    if (verbose)  cmessage('\tScan cas9seqs with %d mismatches against genome', 
                            mismatch)
    vres <- Biostrings::vcountPDict(Biostrings::DNAStringSet(cas9seqs),
                                    bsgenome,
                                    min.mismatch = mismatch,
                                    max.mismatch = mismatch)
    hits <- vres %>%
            data.table::as.data.table() %>%
            extract(, .(n = sum(count)), by = 'index') %>%
            extract2('n')
    
    if (verbose)  report_scan(hits)
    
    return(hits)
}
