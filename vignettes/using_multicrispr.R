## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "    #", 
  cache    = TRUE
)

## ----Flank-----------------------------------------------------------------

targetranges %>% left_flank()
targetranges %>% right_flank()
targetranges %>% slop()
targetranges %<>% flank_fourways()

