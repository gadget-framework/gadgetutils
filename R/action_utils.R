#' @export
g3_action_order <- function(action_name){
  
  action_lookup <- list(
    'time' = 0,
    'printing' = 1,
    'migration' = 2,
    'consumption' = 3,
    'natural_mortality' = 4,
    'growth' = 5,
    'spawning' = 6,
    'maturation' = 7,
    'recruitment' = 8,
    'straying' = 9,
    'likelihood' = 10,
    'reporting' = 11,
    'ageing' = 12
  )
  
  if (!is.character(action_name)){ stop("The 'action_name' argument should be a character string") }
  else{
    if (!casefold(action_name) %in% names(action_lookup)){
      stop("The 'action_name' argument does not match the list of action names")
    }
    return(action_lookup[[match(casefold(action_name), names(action_lookup))]])
  }
}

#' @export
locate_g3_actions <- function(action_list, type, model_object = NULL){
  
  if (!is.list(action_list)) stop("The 'action_list' argument should be a list")
  if (!is.vector(type) || !is.character(type) || length(type) > 2){
    stop("The 'type' parameter should be a character vector with length <= 2")
  }
  
  ## Create a dummy model object is not specified
  if (is.null(model_object)){ model_object <- list(name = NULL) }
  else{
    if (is.character(model_object)){ model_object <- list(name = model_object) }
    else{
      if (!is.list(model_object)){
        stop("The 'model_object' parameter should be a list or character string (of length 1)")
      }
      else{
        if (is.list(model_object) && !("name" %in% names(model_object))){
          stop("The 'model_object' list should contain a 'name' element")
        } 
      }
    }
  }
  
  ## Action types based on the order of calculations:
  ## https://gadget-framework.github.io/gadget2/userguide/chap-order.html
  action_types <- 
    c("time", "initial", "printing", "migration", "consumption", 
      "natural_mortality", "growth", "spawning", "maturation", 
      "recruitment", "straying", "likelihood", "reporting", "ageing")
  
  type_ind <- match(casefold(type), action_types)
  if (any(is.na(type_ind))){
    stop(paste("The following are not valid action types: ", paste(type[which(is.na(type_ind))], collapse = ", ")))
  }
  
  ## Create a look-up table for action_types
  action_lookup <- 
    data.frame(
      ref = action_types,
      code = stringr::str_pad(c(0,0:12), width = 3, side = "left", pad = "0")
    ) 
  
  ## Setup the pattern to look-up
  ## TIME
  if (type == "time"){ action_pattern <- "^000$" }
  else{
    ## OTHER
    action_pattern <- 
      paste("^", 
            paste(action_lookup[type_ind, "code"], collapse = ":"),
            ":", model_object$name, sep = "")
  }
  
  ## Loop over action_list to search for action_pattern 
  out <- list()
  for (i in seq_along(action_list)){
    if (length(action_list[[i]]) == 0) next
    name_ind <- grepl(action_pattern, names(action_list[[i]]))
    if (any(name_ind)){
      tmp <- list(c(pos = as.character(i), 
                    name = names(action_list[[i]])[name_ind], 
                    pattern = action_pattern))
      out <- c(out, tmp)
    }
  }
  
  if (length(out) == 0){ 
    warning(paste0("No action(s) found for the pattern: ", action_pattern))
  }
  else{
    #print(paste0("Number of actions found = ", length(out)))
    out <- as.data.frame(do.call("rbind", out))
    out$pos <- as.numeric(out$pos)
  }
  return(out)
}
