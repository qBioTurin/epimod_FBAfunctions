
#------------------------------------------------------------------------------#
#                       FBA_greatmod class definition                          #
#------------------------------------------------------------------------------#

setClass(
  # Set the name for the class
  "FBA_greatmod",
  # Define the slots
  representation(
    S = "matrix",
    lowbnd = "numeric",
    uppbnd = "numeric",
    obj_coef = "numeric",
    react_id = "character",
    met_id = "character",
    dietName = "character",
    dietValues = "data.frame",
    lowbnd_beforeDiet = "numeric",
    uppbnd_beforeDiet = "numeric"
  )
)

#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

FBA_greatmod <- function(S, ub, lb, obj_fun, react_name=NULL, met_name=NULL) {
  if (missing(S) || missing(ub) || missing(lb) || missing(obj_fun) ) {
    stop("Creating an object of class model needs: S, ub, lb, obj_fun!")
  }
  
  obj <- new("FBA_greatmod", S, ub, lb, obj_fun, react_name, met_name)
  
  stopifnot(validObject(obj))
  
  return(obj)
}

## assign the function as the validity method for the class
validityGreatModClass=function(object)
{
  if(!is.matrix(object@S))
    return("S has to be a matrix with: 1) number of rows equal to number of metabolites, and 2) number of columns equal to the number of reactions.")
  
  ncol = length(object@S[1,])
  nrow = length(object@S[, 1])
  
  if(length(object@obj_coef)!=ncol || length(object@uppbnd)!=ncol || length(object@lowbnd)!=ncol || length(object@react_id)!=ncol){
    return("The obj_fun, ub, and lb have different lengths. The number of rections should be equal to the number of columns of S.")
  }
  
  if(length(unique(object@react_id)) != ncol)
    return("The vector of reactions name is different from the number of columns of S, or multiple reactions have the same name (it is not allowed)!")
  
  if(length(unique(object@met_id))!= nrow)
    return("The vector of metabolites name is different from the number of rows of S, or multiple metabolites have the same name (it is not allowed)!")
  
  return(TRUE)
}

setValidity("FBA_greatmod", validityGreatModClass)

##
#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature("FBA_greatmod"),
          definition = function(.Object, S, ub, lb, obj_fun,
                                react_name=NULL,
                                met_name=NULL) {
            
            if ( (!missing(S)) || (!missing(ub)) || (!missing(lb)) || (!missing(obj_fun)) ) {
              
              .Object@S <- as.matrix(S)
              .Object@lowbnd = lb
              .Object@uppbnd = ub
              .Object@obj_coef = obj_fun
              
              ncol = length(.Object@S[1, ])
              nrow = length(.Object@S[, 1])
              
              if(missing(react_name))
                .Object@react_id = paste0("Reaction", 1:ncol)
              else
                .Object@react_id = react_name
              
              if(missing(met_name))
                .Object@met_id = paste0("Metabolite", 1:nrow)
              else
                .Object@met_id = met_name
              
            }
            return(.Object)
          }
)



# create a method to change the objective function
setGeneric(name="setObjFun",
           def=function(theObject,reaction_name)
           {
             standardGeneric("setObjFun")
           }
)

setMethod(f="setObjFun",
          signature="FBA_greatmod",
          definition=function(theObject, reaction_name)
          {
            if(is.null(theObject@react_id))
              return("No reactions name are present in the class.")
            
            id = which(theObject@react_id %in% reaction_name)
            
            if(length(id) == 0)
              return("The reactions passed in input are not present in the model.")
            
            obj_coef = rep(0,length(theObject@obj_coef))
            obj_coef[id] = 1
            theObject@obj_coef = obj_coef
            
            return(theObject)
          }
)

#' @aliases getExchangesR FBA_greatmod-methods
#' @param theObject A `FBA_greatmod` object
#' @docType methods
#' @rdname FBA_greatmod-methods
#' @return The constraints of the reaction.
#' @export

setGeneric(name="getExchangesR",
           def=function(theObject,reaction.name)
           {
             standardGeneric("getExchangesR")
           }
)

setMethod(f="getExchangesR",
          signature=c("FBA_greatmod"),
          definition=function(theObject)
          {
            
            S = theObject@S
            mat_nonzero <- as.data.frame(which(S != 0, arr.ind = T) ) %>%
              dplyr::group_by(col) %>%
              dplyr::filter(length(col) == 1) %>%
              dplyr::ungroup() %>%
              dplyr::select(col)
            
            if(length(mat_nonzero$col) == 0)
              stop("No exchange reaction was found!")
            else
              ExcR = theObject@react_id[mat_nonzero$col]
            
            return(ExcR)
          }
)

#' @aliases getConstraints FBA_greatmod-methods
#' @param theObject A `FBA_greatmod` object
#' @docType methods
#' @rdname FBA_greatmod-methods
#' @return The constraints of the reaction.
#' @export

setGeneric(name="getConstraints",
           def=function(theObject,reaction.name)
           {
             standardGeneric("getConstraints")
           }
)

setMethod(f="getConstraints",
          signature=c("FBA_greatmod","character"),
          definition=function(theObject,reaction.name)
          {
            index.r = which(theObject@react_id == reaction.name)
            if(length(index.r)==0)
              stop("The reaction does not match any reaction name in the model.")
            
            Constraints = c(theObject@lowbnd[index.r],
                            theObject@uppbnd[index.r])
            
            return(Constraints)
          }
)

#' @aliases setConstraints FBA_greatmod-methods
#' @param theObject A `FBA_greatmod` object
#' @param reaction.name Character of the reaction name.
#' @param newConstraints Numeric vector of length equal to 2 containing the new lower and upper constraints
#' @docType methods
#' @rdname FBA_greatmod-methods
#' @return The FBA_greatmod class with the the constraints updated.
#' @export

setGeneric(name="setConstraints",
           def=function(theObject,reaction.name,newConstraints)
           {
             standardGeneric("setConstraints")
           }
)

setMethod(f="setConstraints",
          signature=c("FBA_greatmod","character","numeric"),
          definition=function(theObject,reaction.name,newConstraints)
          {
            index.r = which(theObject@react_id == reaction.name)
            if(length(index.r)==0)
              stop("The reaction does not match any reaction name in the model.")
            
            if(length(unique(newConstraints)) !=2 )
              stop("The parameter newConstraints must be a vector of two different numeric value (the minimum one will be the new lwbnd, the maximum the uppbnd)  ")
            
            theObject@uppbnd[index.r] = max(newConstraints)
            theObject@lowbnd[index.r] = min(newConstraints)
            
            return(theObject)
          }
)

#' @title FBA_greatmod methods.
#'
#' @name FBA_greatmod-methods
#' @aliases setDiet FBA_greatmod-methods
#' @param theObject A `FBA_greatmod` object
#' @param dietf The diet file or data.frame
#' @docType methods
#' @rdname FBA_greatmod-methods
#' @return The FBA_greatmod class with lower bounds updated with diet values.
#' @export
#'
setGeneric(name="setDiet",
           def=function(theObject,dietf)
           {
             standardGeneric("setDiet")
           }
)

setMethod(f="setDiet",
          signature="FBA_greatmod",
          definition=function(theObject,dietf)
          {
            if(!is.data.frame(dietf)){
              if(file.exists(dietf)){
                theObject@dietName = basename(dietf)
                diet = readr::read_delim(dietf,
                                         delim = "\t", escape_double = FALSE,
                                         trim_ws = TRUE)
                if(length(diet[1,]) != 3)
                  return("The file should have three columns: 1) reactions name, 2) lwbnd and 3) uppbwnd.")
              }
              else
              {
                return("The file does not exist.")
              }
            }
            else if(length(dietf[1,]) == 3)
            {
              diet <- dietf
            }
            else
            {
              return("dietf should be a data.frame with three columns: 1) reactions name, 2) lwbnd and 3) uppbwnd.")
            }
            
            colnames(diet) = c("reactionID","bnd1","bnd2")
            
            theObject@dietValues = diet
            reactionComm = diet$reactionID[diet$reactionID %in% model@react_id]
            
            if(length(reactionComm) == 0)
            {
              stop("No reactions in the diet is shared with the ractions name in the FBA model.")
            }
            else
            {
              bnd_index = match(reactionComm, theObject@react_id )
              theObject@lowbnd_beforeDiet = theObject@lowbnd
              theObject@uppbnd_beforeDiet = theObject@uppbnd
              
              for(b in 1:length(reactionComm))
              {
                bnds = diet[diet$reactionID == reactionComm[b],-1]
                theObject@lowbnd[bnd_index[b]] = min(bnds)
                theObject@uppbnd[bnd_index[b]] = max(bnds)
              }
              
              NotComm = diet$reactionID[! diet$reactionID %in% theObject@react_id]
              if(length(NotComm) > 0 )
                print(paste0( length(NotComm), " reactions of ",length(diet$reactionID)," are not present in the FBA model") )
            }
            
            return(theObject)
          }
)


#' @aliases setDiet.name FBA_greatmod-methods
#' @param theObject A `FBA_greatmod` object
#' @param diet_name Character for the diet name
#' @docType methods
#' @rdname FBA_greatmod-methods
#' @return The FBA_greatmod class with the diet name updated.
#' @export

setGeneric(name="setDiet.name",
           def=function(theObject,diet_name)
           {
             standardGeneric("setDiet.name")
           }
)

setMethod(f="setDiet.name",
          signature=c("FBA_greatmod","character"),
          definition=function(theObject,diet_name)
          {
            theObject@dietName = diet_name
            return(theObject)
          }
)

#' @description Write the fba model in file.
#' @name FBA_greatmod-methods
#' @aliases writeFBAfile FBA_greatmod-methods
#' @param theObject A `FBA_greatmod` object
#' @param fba_fname The file name in which the FBA model will be saved. See *epimod::model.generation*.
#' @docType methods
#' @rdname FBA_greatmod-methods
#' @export
#'
setGeneric(name="writeFBAfile",
           def=function(theObject,fba_fname,...)
           {
             standardGeneric("writeFBAfile")
           }
)

setMethod(f="writeFBAfile",
          signature=c("FBA_greatmod","character"),
          definition=function(theObject,fba_fname=NULL,...)
          {
            if(is.null(fba_fname))
              fba_fname = "fba_file"
            
            print(fba_fname)
            
            ReactionsNames <- unlist(theObject@react_id)
            
            theObject@S -> S
            
            ncol = length(S[1, ])
            nrow = length(S[, 1])
            
            ## Da prevedere che uno possa cambiare b (di default 0 per FBA)
            matrix(0, nrow = ncol, ncol = 1) -> b
            #
            as.matrix(theObject@lowbnd) -> lb
            as.matrix(theObject@uppbnd) -> ub
            c(theObject@obj_coef) -> obj
            
            #### Start writing the file to pass to GLPK solver
            
            print(">> [Writing] Initialization of the file ...")
            
            rb <- cbind(b,b)
            cb <- cbind(lb,ub)
            
            print(">> [Writing] Reactions name ...")
            write(paste(ReactionsNames, collapse = " "), file = fba_fname)
            
            ## Da prevedere che uno possa mettere MIN
            
            print(">> [Writing] Objective function ...")
            write(paste0(nrow," ; ",ncol," ; ","GLP_MAX" ),file = fba_fname, append = T)
            write(paste(obj, collapse = " "),file = fba_fname,append = T )
            
            print(">> [Writing] Costraints ...")
            for(i in 1:nrow) {
              write(paste0("GLP_F"," ; ",paste0(rb[i,],collapse = " ; ") ),file = fba_fname ,append = T)
            }
            
            ## Da prevedere che uno possa cambiare i tipi di bounds
            for(j in 1:ncol){
              write(paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")),file = fba_fname ,append = T)
            }
            print(">> [Writing] Stoichiometric matrix ...")
            for(i in 1:nrow) {
              for(j in 1:ncol){
                write(paste0(i," ; ", j, " ; ", S[i,j]),file = fba_fname ,append = T)
              }
            }
            print(">> [Writing] End.")
          }
)
