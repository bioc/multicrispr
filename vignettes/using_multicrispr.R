## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "    #", 
  cache = TRUE
)

## ----score_rs1-------------------------------------------------------------
cas9ranges$score_rs1 <- score_cas9ranges(cas9ranges, ruleset = 1)

## ----score_rs2-------------------------------------------------------------
python_available <- assertive.reflection::is_unix() & 
                    reticulate::py_module_available('azimuth')
if (python_available){
  cas9ranges$score_rs2 <- score_cas9s(
                            cas9ranges, ruleset = 2, condaenv = 'default')
  
  message('\tFilter for scores > 0.4')
  cas9ranges %<>% extract(cas9ranges$score > 0.4)
  cmessage('\t\t%d cas9 seqs across %d ranges', 
            length(unique(cas9ranges$seqs)), 
            length(cas9ranges))
}

