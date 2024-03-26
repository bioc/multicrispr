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
# BASE
#===========

    #-------
    # engine
    #-------

        safe_deparse <- function(expr, ...){
            paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
        }
        
        
        get_name_in_parent <- function(x, escape_percent = TRUE){
            xname <- safe_deparse(  do.call(  substitute, 
                                              list(substitute(x), parent.frame()) )  )
            if(escape_percent)  xname <- gsub("%", "%%", xname)
            xname
        }

        
        cause <- function(x){
            y <- attr(x, "cause")
            if(is.null(y))  return(noquote(character(length(x))))
            y
        }

                
        `cause<-` <- function(x, value){
            # Can't use is_scalar here due to dependency on this
            if(length(value) != 1 && length(value) != length(x)){
                stop(  sprintf( "The length of value should be 1 or the length of x (%d) but is %d.", 
                                length(x),
                                length(value) ) )
            }
            attr(x, "cause") <- noquote(as.character(value))
            x
        }
        
        
        false <- function(...){
            msg <- if(nargs() > 0L) sprintf(...) else ""
            x <- FALSE
            cause(x) <- msg[1]
            class(x) <- c("scalar_with_cause", "logical")
            x
        }


        set_cause <- function(x, false_value, missing_value = "missing"){
            if(!anyNA(x) && all(x, na.rm = TRUE)) return(x)  # fast version of all(!is.na(x) & x)
            is_na_x <- is.na(x)
            len_x <- length(x)
            # TRUES
                cause_value <- character(len_x)
            # NAS
                if(length(missing_value) == 1){
                    cause_value[is_na_x] <- missing_value
                } else {
                    missing_value <- rep_len(missing_value, len_x)
                    cause_value[is_na_x] <- missing_value[is_na_x]
                }
            # FALSES
                false_index <- !(x | is_na_x) # more efficient to calc than !x & !is_na_x
                if(length(false_value) == 1){
                    cause_value[false_index] <- false_value
                } else {
                    false_value <- rep_len(false_value, len_x)
                    cause_value[false_index] <- false_value[false_index]
                }
            cause(x) <- cause_value
            class(x) <- c("vector_with_cause", "logical")
            x
        }

        
        to_names <- function(x){                                 # special handling for double, complex only
            if(is.double(x) && is.vector(x)){
                ifelse(is.na(x), NA_real_, sprintf("%.17g", x))  # is.vector prevents matching to POSIXct
            } else if(is.complex(x)){
                ifelse(is.na(x), NA_complex_, sprintf("%.17g+%.17gi", Re(x), Im(x)))
            } else {
                as.character(x)
            }
        }
        
        
        bapply <- function(x, predicate, ...){
            vapply(x, predicate, logical(1L), ..., USE.NAMES = TRUE)
        }
        
        
        call_and_name <- function(fn, x, ...){
            y <- fn(x, ...)
            dim(y) <- dim(x)
            names(y) <- to_names(x)
            y
        }
        
        
        give_feedback <- function(handler_type, msg, predicate_name){
            handler <- match.fun(handler_type)
            ass_condition <- switch(
                handler_type,
                stop = assertionError,
                warning = assertionWarning,
                message = assertionMessage)
            # Throw error/warning/message
            caller <- if(sys.nframe() >= 3){  sys.call(-3)
            } else {                          NULL         }
            
            # UTF-8 characters do not display correctly under Windows for some 
            # LC_CTYPE locale values, but there isn't much assertive can do about that.
            # https://stackoverflow.com/q/32696241/134830
            handler(ass_condition(paste(predicate_name, msg, sep = " : "), caller, predicate_name))
        }
        
        
        assert_engine <- function(
            predicate, 
                  ...,
                  msg = "The assertion failed.",
                 what = c("all", "any"),
            na_ignore = FALSE,
             severity = c("stop", "warning", "message", "none")
        ){
            handler_type <- match.arg(severity)
            dots <- list(...)
            return_value <- if(length(dots) > 0) dots[[1]] else NULL
            if(handler_type == "none")  return(invisible(return_value))
            what <- match.fun(match.arg(what))
            predicate_name <- get_name_in_parent(predicate)
            
            ok <- predicate(...)
            if(inherits(ok, "scalar_with_cause")){
                if(!isTRUE(ok)){
                    if(missing(msg))   msg <- cause(ok)
                    give_feedback(handler_type, msg, predicate_name)
                }
            } else { # inherits(ok, "vector_with_cause")
                really_ok <- if(na_ignore){ ok |  is.na(ok)   # ok can be TRUE or NA; FALSE is bad
                             } else {       ok & !is.na(ok) } # ok can be TRUE; FALSE or NA is bad
                if(!what(really_ok)){
                    # Append first few failure values and positions to the error message.
                    msg <- paste(enc2utf8(msg), print_and_capture(ok), sep = "\n")
                    give_feedback(handler_type, msg, predicate_name)
                }
            }
            invisible(return_value)
        }
        
        
        is2 <- function(x, class, .xname = get_name_in_parent(x)){    
            # Can't use is_empty in next line because that function calls this one.
            if(length(class) == 0L) stop("You must provide a class.")
            if(length(class) > 1L)  return(  set_cause( bapply(class, function(cl) is2(x, cl, "")),
                                                        sprintf("%s is not '%s'", 
                                                                type_description(x), class))  )
            ok <- tryCatch(  {   is.class <- match.fun(paste0("is.", class))
                                 is.class(x)
                             },
                             error = function(e) is(x, class)  )
            if(!ok) return( false( "%s is not of class '%s'; it has %s.", 
                                   .xname, class, type_description(x)  )  )
            TRUE
        }        
        
        
        coerce_to <- function(x, target_class, .xname = get_name_in_parent(x)){
            # Can't use is_empty in next line because that function calls this one.
            if(length(target_class) == 0L)   stop("You must provide a class.")
            if(!is.character(target_class))  stop("target_class should be a character vector.")
            for(this_class in target_class){
                if(!is2(x, this_class)){
                    warning(sprintf("Coercing %s to class %s.", .xname, sQuote(this_class)),
                            call. = FALSE)  }
                tryCatch(  
                    {   as.this_class <- match.fun(paste0("as.", this_class))
                        return(as.this_class(x))
                    },
                    error = function(e){
                        # as.this_class doesn't exist; try as(, "this_class") instead
                        tryCatch(   return(as(x, this_class)),
                                    error = function(e){
                                        # Can't coerce to this class; warn and move to next class
                                        warning( sprintf("%s cannot be coerced to type %s.", 
                                                         .xname, sQuote(this_class)), 
                                                 call. = FALSE )
                                    } ) } )
            }
            # Nothing worked; throw an error
            stop(sprintf("%s cannot be coerced to any of these types: %s.", 
                         .xname, toString(sQuote(target_class))))
        }

                
        type_description <- function(x){
          if(is.array(x)){           sprintf(sprintf("class '%s %s'", class(x[FALSE]), toString(class(x))))
          } else if(is.function(x)){ sprintf(sprintf("class '%s %s'", typeof(x), toString(class(x))))
          } else if(isS4(x)){        sprintf(sprintf("S4 class '%s'", toString(class(x))))
          } else {                   sprintf("class '%s'", toString(class(x)))
          }
        }


        strip_attributes <- function(x){
            attributes(x) <- NULL
            x
        }
        
        
        use_first <- function(x, indexer = c("[[", "["), .xname = get_name_in_parent(x)){
            len_x <- length(x)
            # Can't use assert_is_non_empty, is_scalar in next lines because those 
            # functions calls this one.
            if(len_x == 0L)  stop(sprintf("%s has length 0.", .xname))
            if(len_x == 1L)  return(x)
            indexer <- match.fun(match.arg(indexer))
            x1 <- indexer(x, 1L)
            warning( sprintf("Only the first value of %s (= %s) will be used.", 
                             .xname, as.character(x1)),
                     call. = FALSE )
            x1
        }
        
        
    #-------------
    # true / false
    #-------------

        
        is_identical_to_true <- function(
            x, allow_attributes = FALSE, .xname = get_name_in_parent(x)
        ){
            if(allow_attributes)  x <- strip_attributes(x)
            if(!identical(TRUE, x)){
                msg <- gettextf(  "%s is not identical to TRUE; its value is %s.", 
                                  .xname, 
                                   safe_deparse(x),
                                   domain = "R-assertive.base"  )
                return(false(msg))
            }
            TRUE
        }
        
        
        is_true <- function(x, .xname = get_name_in_parent(x)){
            x <- coerce_to(x, "logical", .xname)
            call_and_name(  function(x){  is_na_x <- is.na(x)
                                          ok <- x & !is_na_x
                                          set_cause(ok, ifelse(is_na_x, "missing", "false"))  }, 
                            x  )
        }
        
        
        is_false <- function(x, .xname = get_name_in_parent(x)){
            x <- coerce_to(x, "logical", .xname)
            call_and_name(  function(x){   is_na_x <- is.na(x)
                                           ok <- !(x | is_na_x)  # same as !x & !is_na_x
                                           set_cause(ok, ifelse(is_na_x, "missing", "true"))  }, 
                            x  )
        }
        

        assert_all_are_true <- function(x, severity = getOption("assertive.severity", "stop")){                                                     
            msg <- gettextf(  "The values of %s are not all TRUE.", 
                               get_name_in_parent(x), 
                               domain = "R-assertive.base"  )
            assert_engine(  is_true, 
                                  x, 
                                msg = msg, 
                             .xname = get_name_in_parent(x), 
                           severity = severity )
        }
        

        assert_all_are_false <- function(x, severity = getOption("assertive.severity", "stop")){                                                     
            msg <- gettextf(  "The values of %s are not all FALSE.", 
                              get_name_in_parent(x), 
                              domain = "R-assertive.base"  )
            assert_engine(  is_false, 
                                   x, 
                                 msg = msg, 
                              .xname = get_name_in_parent(x), 
                            severity = severity  )
        }

        
    #----------
    # identical
    #----------

        are_identical <- function(
                           x, 
                           y, 
            allow_attributes = FALSE, 
                      .xname = get_name_in_parent(x),
                      .yname = get_name_in_parent(y)
        ){  
            if(allow_attributes){
                x <- strip_attributes(x)
                y <- strip_attributes(y)
            }
            if(!identical(x, y))  return( false( gettext("%s and %s are not identical."), 
                                                 .xname, .yname ) )
            TRUE
        }
        
        
        assert_are_identical <- function(
            x, y, allow_attributes = FALSE, severity = getOption("assertive.severity", "stop")
        ){
            assert_engine(  are_identical,
                                        x, 
                                        y = y,
                                   .xname = get_name_in_parent(x),
                                   .yname = get_name_in_parent(y),
                                 severity = severity  )
        }
        
        
        
#===========
# PROPERTIES
#===========


    #-----------
    # duplicates
    #-----------


        has_duplicates <- function(x, .xname = get_name_in_parent(x)){
            if(!anyDuplicated(x))  return(false(gettext("%s has no duplicates."), .xname))
            TRUE
        }
        
        
        has_no_duplicates <- function(x, .xname = get_name_in_parent(x)){
            if(anyDuplicated(x)){
                dupe_indicies <- which(duplicated(x))
                return( false( ngettext( length(dupe_indicies),
                                     "%s has a duplicate at position %s.",
                                     "%s has duplicates at positions %s." ), 
                              .xname, 
                               toString(dupe_indicies, width = 100)  )  )
            }
            TRUE
        }

        
        assert_has_duplicates <- function(x, severity = getOption("assertive.severity", "stop")){                                                                
            assert_engine(  has_duplicates, 
                            x, 
                            .xname = get_name_in_parent(x), 
                            severity = severity )
        }
        
        
        assert_has_no_duplicates <- function(x, severity = getOption("assertive.severity", "stop")){
            assert_engine(  has_no_duplicates,
                            x,
                           .xname = get_name_in_parent(x),
                            severity = severity )
        }
        
        
    #-------
    # length
    #-------
        
        are_same_length <- function(
            x, y, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y)
        ){
            len_x <- length(x)
            len_y <- length(y)
            if(len_x != len_y){
                return( false(  gettext("%s has length %d but %s has length %d."),
                               .xname,
                                len_x,
                               .yname,
                                len_y  ) )
            }
            TRUE
        }
        
        assert_are_same_length <- function(
            x, y,  severity = getOption("assertive.severity", "stop")){
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
        
        is_of_length <- function(x, n, .xname = get_name_in_parent(x)){
            n <- use_first(n)
            check_n(n)
            length_x <- length(x)
            if(length_x != n)  return(false("%s has length %d, not %d.", .xname, length_x, n))
            TRUE
        }

        
        DIM <- function(x){
            dim_x <- dim(x)
            if(is.null(dim_x)) length(x) else dim_x
        }
        
        
        n_elements <- function(x){
            if(is.recursive(x)){    sum(vapply(x, n_elements, integer(1)))
            } else {                as.integer(prod(DIM(x)))
            }  
        }
        
                
        has_elements <- function(x, n, .xname = get_name_in_parent(x)){
            n <- use_first(n)
            check_n(n)
            n_elements_x <- n_elements(x)
            if(n_elements_x != n){
                return(  false(  ngettext( n_elements_x, 
                                          "%s has %d element, not %d.", 
                                          "%s has %d elements, not %d."),
                                .xname, 
                                 n_elements_x,
                                 n ) )
            }
            TRUE
        }
        
        
        get_metric <- function(metric){
            switch(  metric,
                     length = is_of_length,
                   elements = has_elements,
                   stop("Bug in assertive; the metric", metric, "is not valid.", domain = NA) )
        }
    
        
        is_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x)){
            metric <- match.arg(metric)
            metric_fn <- get_metric(metric)
            metric_fn(x, 1L, .xname)
        }     
        
        
        is_non_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x)){
            metric <- match.arg(metric)
            metric_fn <- get_metric(metric)
            if(metric_fn(x, 1)){
                msg <- switch(  metric,
                                length = gettext("%s has length 1."),
                              elements = gettext("%s has 1 element.") )
                return(false(msg, .xname))
            }
            TRUE
        }
        
        
        assert_is_scalar <- function(
            x, metric = c("length", "elements"), severity = getOption("assertive.severity", "stop")
        ){                                        
            metric <- match.arg(metric)
            assert_engine(  is_scalar, 
                                    x, 
                               metric = metric, 
                               .xname = get_name_in_parent(x),
                             severity = severity  )
        }
        
        
        assert_is_non_scalar <- function(
            x, metric = c("length", "elements"), severity = getOption("assertive.severity", "stop")
        ){                            
            metric <- match.arg(metric)
            assert_engine(  is_non_scalar,
                                        x,
                                   metric = metric,
                                   .xname = get_name_in_parent(x),
                                 severity = severity  )
        }
        
        
        is_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x)){  
            metric <- match.arg(metric)
            metric_fn <- get_metric(metric)
            metric_fn(x, 0L, .xname)
        }        
        
    #-------
    # names
    #-------
        
        has_names <- function(x, .xname = get_name_in_parent(x)){
            namesx <- names(x)
            if(    is.null(namesx))   return(false("The names of %s are NULL.",      .xname))
            if(!any(nzchar(namesx)))  return(false("The names of %s are all empty.", .xname))
            TRUE
        } 


        assert_has_names <- function(x, severity = getOption("assertive.severity", "stop")){                                                            
            assert_engine(  has_names, 
                                    x, 
                               .xname = get_name_in_parent(x),
                             severity = severity  )
        }
        
                        
#======
# TYPES
#======

    #--------
    # logical
    #--------
        
        is_logical <- function(x, .xname = get_name_in_parent(x)){
            is2(x, "logical", .xname)
        }       
        
        
        is_a_bool <- function(x, .xname = get_name_in_parent(x)){
            if(!(ok <- is_logical(x, .xname)))          return(ok)
            if(!(ok <- is_scalar( x, .xname = .xname))) return(ok)
            TRUE
        }
        
        
        assert_is_a_bool <- function(x, severity = getOption("assertive.severity", "stop")){      
            assert_engine(  is_a_bool, 
                                    x, 
                               .xname = get_name_in_parent(x), 
                             severity = severity  )
        }
        
        
    #-------
    # number
    #-------
        
        
        is_numeric <- function(x, .xname = get_name_in_parent(x)){
            is2(x, "numeric", .xname)
        }
        
        
        is_a_number <- function(x, .xname = get_name_in_parent(x)){
            if(!(ok <- is_numeric(x, .xname)))           return(ok)
            if(!(ok <- is_scalar( x, .xname = .xname)))  return(ok)
            TRUE
        } 
        
        
        assert_is_a_number <- function(x, severity = getOption("assertive.severity", "stop")){                                                          
            assert_engine(  is_a_number, 
                                      x, 
                                 .xname = get_name_in_parent(x), 
                               severity = severity  )
        }
    
        
        assert_is_numeric <- function(x, severity = getOption("assertive.severity", "stop")){                                                         
            assert_engine(  is_numeric, 
                                     x, 
                                .xname = get_name_in_parent(x),
                              severity = severity  )
        }

                
    #-------------------
    # character / string
    #-------------------

                
        is_character <- function(x, .xname = get_name_in_parent(x)){
            is2(x, "character", .xname)
        }
        
        
        is_a_string <- function(x, .xname = get_name_in_parent(x)){
            if(!(ok <- is_character(x, .xname)))           return(ok)
            if(!(ok <- is_scalar(   x, .xname = .xname)))  return(ok)
            TRUE
        }
        
        
        assert_is_character <- function(x, severity = getOption("assertive.severity", "stop")){                                                         
            assert_engine(  is_character, 
                                       x, 
                                  .xname = get_name_in_parent(x),
                                severity = severity  )
        }

                
        assert_is_a_string <- function(x, severity = getOption("assertive.severity", "stop")){                                                         
            assert_engine(  is_a_string, 
                                      x, 
                                 .xname = get_name_in_parent(x), 
                               severity = severity  )
        }
        
        
    #------
    # class
    #------
    
        assert_is_all_of <- function(
            x, classes, severity = getOption("assertive.severity", "stop")
        ){  
            msg <- gettextf(  "%s is not in all of the classes %s.", 
                               get_name_in_parent(x), 
                               toString(sQuote(classes)) )
            assert_engine(  is2, 
                              x, 
                          class = classes, 
                            msg = msg, 
                       severity = severity  )
        }

        
        assert_is_any_of <- function(
            x, classes, severity = getOption("assertive.severity", "stop")){  
            msg <- gettextf(  "%s is not in any of the classes %s.", 
                               get_name_in_parent(x), 
                               toString(sQuote(classes)) )
            assert_engine(  is2, 
                              x, 
                          class = classes, 
                            msg = msg, 
                           what = "any",
                       severity = severity  )
        }
        
        
#=========
# STRINGS
#=========

    #------
    # regex
    #------

        is_matching_regex <- function(
            x, pattern, opts_regex = NULL, .xname = get_name_in_parent(x)
        ){
            x <- coerce_to(x, "character", .xname)
            call_and_name(
                function(x){
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
            msg <- sprintf( "%s does not match %s", .xname, .fixedname )
            assert_engine(  is_matching_regex, 
                                            x,
                                      pattern,
                                   opts_regex = opts_regex,
                                       .xname = .xname,
                                          msg = msg, 
                                         what = 'all',
                                    na_ignore = na_ignore,
                                     severity = severity
            )
        }
        

#===========
# REFLECTION
#===========


    #-----------
    # windows
    #-----------

        not_this_os <- function(os){
            false(  
                gettext("The operating system is not %s. R reports it as: Sys.info()['sysname'] = %s, .Platform$OS = %s."), 
                os, 
                Sys.info()["sysname"],
                .Platform$OS
            )
        }
    
    
        is_windows <- function(){
            if(.Platform$OS.type != "windows")  return(not_this_os("Windows"))
            TRUE
        }
    

