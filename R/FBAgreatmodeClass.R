
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
    dietName = "character",
    dietValues = "data.frame",
    lowbnd_beforeDiet = "numeric"
  )
)

#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

FBA_greatmod <- function(S, ub, lb, obj_fun, react_name=NULL) {
  if (missing(S) || missing(ub) || missing(lb) || missing(obj_fun) ) {
    stop("Creating an object of class model needs: S, ub, lb, obj_fun!")
  }

  obj <- new("FBA_greatmod", S, ub, lb, obj_fun, react_name)

  stopifnot(validObject(obj))

  return(obj)
}

## assign the function as the validity method for the class
validityGreatModClass=function(object)
{
  if(!is.matrix(object@S))
    return("S has to be a matrix with: 1) number of rows equal to number of metabolites, and 2) number of columns equal to the number of reactions.")

  ncol = length(object@S[1,])

  if(length(object@obj_coef)!=ncol || length(object@uppbnd)!=ncol || length(object@lowbnd)!=ncol || length(object@react_id)!=ncol){
    return("The obj_fun, ub, and lb have different lengths. The number of rections should be equal to the number of columns of S.")
  }

  if(length(object@react_id)!=ncol)
    return("The vector of reactions id is different from the number of columns of S!")

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
                                react_name=NULL) {

            if ( (!missing(S)) || (!missing(ub)) || (!missing(lb)) || (!missing(obj_fun)) ) {

              .Object@S <- as.matrix(S)
              .Object@lowbnd = lb
              .Object@uppbnd = ub
              .Object@obj_coef = obj_fun

              ncol = length(.Object@S[1,])

              if(missing(react_name))
                .Object@react_id = paste0("Flux",1:ncol)
              else
                .Object@react_id = react_name
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
          definition=function(theObject,reaction_name)
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
                if(length(diet[1,]) != 2)
                  return("The file should have two columns: 1) reactions name and 2) flux values.")
              }
              else
              {
                return("The file does not exist.")
              }
            }
            else if(length(dietf[1,]) == 2)
            {
              diet <- dietf
            }
            else
            {
              return("dietf should be a data.frame with two columns: 1) reactions name and 2) flux values.")
            }

            colnames(diet) = c("reactionID","fluxValue")

            if(length(which(diet$fluxValue > 0)))
            {
              diet$fluxValue[which(diet$fluxValue > 0)] = -diet$fluxValue[which(diet$fluxValue > 0)]
            }

            theObject@dietValues = diet
            reactionComm = diet$reactionID[diet$reactionID %in% model@react_id]

            if(length(reactionComm) == 0)
            {
              stop("No reactions in the diet is shared with the ractions name in the FBA model.")
            }
            else
            {
              lowbnd_index = match(reactionComm, theObject@react_id )
              theObject@lowbnd_beforeDiet = theObject@lowbnd

              for(b in 1:length(reactionComm))
              {
                theObject@lowbnd[lowbnd_index[b]] = diet$fluxValue[diet$reactionID == reactionComm[b]]
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

            ncol=length(S[1,])
            nrow=length(S[,1])

            ## Da prevedere che uno possa cambiare b (di default 0 per FBA)
            matrix(0,nrow = ncol, ncol = 1) -> b
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
