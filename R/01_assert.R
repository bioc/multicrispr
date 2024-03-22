#----------------------------------------------------------------------------------------------
#
#                      `assertive` functions - Richard Cotton
#
#     `assertive` was a suite of CRAN packages by Richard Cotton.
#      Its usage was beautifully documented in the O Reilly book `Testing R` by Richard Cotton.
#      The suite is still fully available at bitbucket/richierocks.
#      But, sadly, the CRAN packages are being deprecated.
#      Multicrispr, as BioC package, requires dependencies to be available on CRAN (or BioC).
#      Therefore, in response, used assertive functionality is now copied into this file.
#      Explicit dependency on assertive is being phased out.
#      Gratefulness towards Richard Cotton for his amazing functionality remains : )
#
#----------------------------------------------------------------------------------------------


#===========
# PROPERTIES
#===========


    #-----------
    # duplicates
    #-----------


        has_duplicates <- function(x, .xname = get_name_in_parent(x))
        {
            if(!anyDuplicated(x)) 
            {
                return(false(gettext("%s has no duplicates."), .xname))
            }
            TRUE
        }
        
        
        has_no_duplicates <- function(x, .xname = get_name_in_parent(x))
        {
            if(anyDuplicated(x)) 
            {
                dupe_indicies <- which(duplicated(x))
                return(
                    false(
                        ngettext(
                            length(dupe_indicies),
                            "%s has a duplicate at position %s.",
                            "%s has duplicates at positions %s."
                        ), 
                        .xname, 
                        toString(dupe_indicies, width = 100)
                    )
                )
            }
            TRUE
        }

        
        assert_has_duplicates <- function(x, 
                                          severity = getOption("assertive.severity", "stop"))
        {                                                                
            assert_engine(
                has_duplicates, 
                x, 
                .xname = get_name_in_parent(x),
                severity = severity
            )
        }
        assert_has_no_duplicates <- function(x, severity = getOption("assertive.severity", "stop"))
        {                                                             
            assert_engine(
                has_no_duplicates, 
                x, 
                .xname = get_name_in_parent(x),
                severity = severity
            )
        }
        
    #-------
    # length
    #-------
        
        are_same_length <- function(x, y, .xname = get_name_in_parent(x),
                                    .yname = get_name_in_parent(y))
        {
            len_x <- length(x)
            len_y <- length(y)
            if(len_x != len_y)
            {
                return(
                    false(
                        gettext("%s has length %d but %s has length %d."),
                        .xname,
                        len_x,
                        .yname,
                        len_y
                    )
                )
            }
            TRUE
        }
        
        assert_are_same_length <- function(x, y, 
                                           severity = getOption("assertive.severity", "stop"))
        {
            assert_engine(
                are_same_length,
                x, 
                y = y,
                .xname = get_name_in_parent(x),
                .yname = get_name_in_parent(y),
                severity = severity
            )
        }
        

        check_n <- function(n)
        {
            if(any(n < 0 | n != round(n)))
            {
                stop("n should be a non-negative integer vector.")
            }
        }
        
        is_of_length <- function(x, n, .xname = get_name_in_parent(x))
        {
            n <- use_first(n)
            check_n(n)
            length_x <- length(x)
            if(length_x != n)
            {
                return(false("%s has length %d, not %d.", .xname, length_x, n))
            }
            TRUE
        }

        
        DIM <- function(x)
        {
            dim_x <- dim(x)
            if(is.null(dim_x)) length(x) else dim_x
        }
        
        
        n_elements <- function(x)
        {
            if(is.recursive(x))
            {
                sum(vapply(x, n_elements, integer(1)))
            } else
            {
                as.integer(prod(DIM(x)))
            }  
        }
        
                
        has_elements <- function(x, n, .xname = get_name_in_parent(x))
        {
            n <- use_first(n)
            check_n(n)
            n_elements_x <- n_elements(x)
            if(n_elements_x != n)
            {
                return(
                    false(
                        ngettext(
                            n_elements_x, 
                            "%s has %d element, not %d.", 
                            "%s has %d elements, not %d."
                        ),
                        .xname, 
                        n_elements_x,
                        n
                    )
                )
            }
            TRUE
        }
        
        get_metric <- function(metric)
        {
            switch(
                metric,
                length   = is_of_length,
                elements = has_elements,
                stop("Bug in assertive; the metric", metric, "is not valid.", domain = NA)
            )
        }
    
        
        is_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
        {
            metric <- match.arg(metric)
            metric_fn <- get_metric(metric)
            metric_fn(x, 1L, .xname)
        }     
        
        
        is_non_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
        {
            metric <- match.arg(metric)
            metric_fn <- get_metric(metric)
            if(metric_fn(x, 1)) 
            {
                msg <- switch(
                    metric,
                    length = gettext("%s has length 1."),
                    elements = gettext("%s has 1 element.")
                )
                return(false(msg, .xname))
            }
            TRUE
        }
        
        assert_is_scalar <- function(x, metric = c("length", "elements"), 
                                     severity = getOption("assertive.severity", "stop"))
        {                                        
            metric <- match.arg(metric)
            assert_engine(
                is_scalar, 
                x, 
                metric = metric, 
                .xname = get_name_in_parent(x),
                severity = severity
            )
        }        
        assert_is_non_scalar <- function(x, metric = c("length", "elements"), 
                                         severity = getOption("assertive.severity", "stop"))
        {                            
            metric <- match.arg(metric)                                 
            assert_engine(
                is_non_scalar, 
                x, 
                metric = metric, 
                .xname = get_name_in_parent(x),
                severity = severity
            )   
        }
        
        
        is_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
        {  
            metric <- match.arg(metric)
            metric_fn <- get_metric(metric)
            metric_fn(x, 0L, .xname)
        }        
        
    #-------
    # names
    #-------
        
        has_names <- function(x, .xname = get_name_in_parent(x))
        {
            namesx <- names(x)
            if(is.null(namesx)) 
            {
                return(false("The names of %s are NULL.", .xname))
            }
            if(!any(nzchar(namesx))) 
            {
                return(false("The names of %s are all empty.", .xname))
            }
            TRUE
        } 


        assert_has_names <- function(x, severity = getOption("assertive.severity", "stop"))
        {                                                            
            assert_engine(
                has_names, 
                x, 
                .xname = get_name_in_parent(x),
                severity = severity
            )
        }
        
                        
#======
# TYPES
#======

    #--------
    # logical
    #--------
        
        is_logical <- function(x, .xname = get_name_in_parent(x))
        {
            is2(x, "logical", .xname)
        }       
        
        
        is_a_bool <- function(x, .xname = get_name_in_parent(x))
        {
            if(!(ok <- is_logical(x, .xname))) return(ok)
            if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
            TRUE
        }
        
        
        assert_is_a_bool <- function(x, severity = getOption("assertive.severity", "stop"))
        {      
            assert_engine(
                is_a_bool, 
                x, 
                .xname = get_name_in_parent(x), 
                severity = severity
            )
        }
        
    #-------
    # number
    #-------
        
        is_numeric <- function(x, .xname = get_name_in_parent(x))
        {
            is2(x, "numeric", .xname)
        }
        
        
        is_a_number <- function(x, .xname = get_name_in_parent(x))
        {
            if(!(ok <- is_numeric(x, .xname))) return(ok)
            if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
            TRUE
        } 
        
        
        assert_is_a_number <- function(x, severity = getOption("assertive.severity", "stop"))
        {                                                          
            assert_engine(
                is_a_number, 
                x, 
                .xname = get_name_in_parent(x), 
                severity = severity
            )
        }
    
        assert_is_numeric <- function(x, 
                                      severity = getOption("assertive.severity", "stop"))
        {                                                         
            assert_engine(
                is_numeric, 
                x, 
                .xname = get_name_in_parent(x),
                severity = severity
            )
        }
        
        #-------------------
    # character / string
    #-------------------
        
        is_character <- function(x, .xname = get_name_in_parent(x))
        {
            is2(x, "character", .xname)
        }
        
        
        is_a_string <- function(x, .xname = get_name_in_parent(x))
        {
            if(!(ok <- is_character(x, .xname))) return(ok)
            if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
            TRUE
        }
        
        assert_is_character <- function(x, 
                                        severity = getOption("assertive.severity", "stop"))
        {                                                         
            assert_engine(
                is_character, 
                x, 
                .xname = get_name_in_parent(x),
                severity = severity
            )
        }

                
        assert_is_a_string <- function(x, 
                                       severity = getOption("assertive.severity", "stop"))
        {                                                         
            assert_engine(
                is_a_string, 
                x, 
                .xname = get_name_in_parent(x), 
                severity = severity
            )
        }
        
        
    #------
    # class
    #------
    
        assert_is_all_of <- function(x, classes, severity = getOption("assertive.severity", "stop"))
        {  
            msg <- gettextf(
                "%s is not in all of the classes %s.", 
                get_name_in_parent(x), 
                toString(sQuote(classes))
            )
            assert_engine(
                is2, 
                x, 
                class = classes, 
                msg = msg, 
                severity = severity
            )
        }    

        
        assert_is_any_of <- function(x, classes, severity = getOption("assertive.severity", "stop"))
        {  
            msg <- gettextf(
                "%s is not in any of the classes %s.", 
                get_name_in_parent(x), 
                toString(sQuote(classes))
            )
            assert_engine(
                is2, 
                x, 
                class = classes, 
                msg = msg, 
                what = "any",
                severity = severity
            )
        }
        
        
#=========
# STRINGS
#=========

    #------
    # regex
    #------

        is_matching_regex <- function(x, pattern, opts_regex = NULL, .xname = get_name_in_parent(x))
        {
            x <- coerce_to(x, "character", .xname)
            call_and_name(
                function(x)
                {
                    ok <- stringi::stri_detect_regex(x, pattern, opts_regex = opts_regex)
                    set_cause(ok, gettextf("does not match '%s'", pattern))
                },
                x
            )
        }
        
        
        assert_all_are_matching_regex <- function(
            x, pattern, opts_regex = NULL, na_ignore = FALSE, severity = getOption("assertive.severity", "stop")
        ){
            .xname <- get_name_in_parent(x)                                             
            .fixedname <- get_name_in_parent(pattern)
            msg <- sprintf(
                "%s does not match %s", 
                .xname, .fixedname
            )
            assert_engine(
                is_matching_regex, 
                x,
                pattern,
                opts_regex = opts_regex,
                .xname     = .xname,
                msg        = msg, 
                what       = 'all',
                na_ignore  = na_ignore,
                severity   = severity
            )
        }
        

#===========
# REFLECTION
#===========


    #-----------
    # windows
    #-----------

        not_this_os <- function(os)
        {
            false(
                gettext(
                    "The operating system is not %s. R reports it as: Sys.info()['sysname'] = %s, .Platform$OS = %s."
                ), 
                os, 
                Sys.info()["sysname"],
                .Platform$OS
            )
        }
    
    
        is_windows <- function()
        {
            if(.Platform$OS.type != "windows")
            {
                return(not_this_os("Windows"))
            }
            TRUE
        }
    

