#' An S4 Class for Flux Balance Analysis Models
#'
#' The `FBA_greatmod` class stores all relevant information to perform Flux Balance Analysis,
#' including stoichiometric matrix, reaction bounds, objective coefficients, and annotations.
#' @rdname FBA_greatmod-class
#' @slot S Numeric matrix. Stoichiometric matrix (metabolites × reactions).
#' @slot lowbnd Numeric vector. Lower bounds of reactions.
#' @slot uppbnd Numeric vector. Upper bounds of reactions.
#' @slot obj_coef Numeric vector. Objective function coefficients.
#' @slot react_id Character vector. Reaction identifiers.
#' @slot met_id Character vector. Metabolite identifiers.
#' @slot dietName Character string. Name of the diet used in simulations.
#' @slot dietValues Data frame. Contains dietary constraint values.
#' @slot lowbnd_beforeDiet Numeric vector. Lower bounds before dietary constraints.
#' @slot uppbnd_beforeDiet Numeric vector. Upper bounds before dietary constraints.
#' @slot bioMax Numeric. Maximum biomass yield (optional).
#' @slot bioMean Numeric. Mean biomass yield (optional).
#' @slot bioMin Numeric. Minimum biomass yield (optional).
#' @slot gene_assoc Numeric vector. Gene associations per reaction.
#' @slot pFBAFlag Numeric. Flag to indicate parsimonious FBA configuration.
#'
#' @export

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
    uppbnd_beforeDiet = "numeric",
    #biomass values
    bioMax = "numeric",
    bioMean = "numeric",
    bioMin = "numeric",
    gene_assoc = "numeric",
    pFBAFlag = "numeric"
  )
)

#' Create an FBA_greatmod Object
#'
#' Constructs a new object of the `FBA_greatmod` S4 class using model inputs.
#'
#' @inheritParams FBA_greatmod-class
#'
#' @return An object of class `FBA_greatmod`
#' @export

FBA_greatmod <- function(S, ub, lb, obj_fun, react_name=NULL, met_name=NULL, bioMax = -1, bioMean = -1, bioMin = -1, gene_assoc = NULL, pFBAFlag = -1) {
  if (missing(S) || missing(ub) || missing(lb) || missing(obj_fun) ) {
    stop("Creating an object of class model needs: S, ub, lb, obj_fun!")
  }

  obj <- new("FBA_greatmod",
             S, ub, lb, obj_fun,
             react_name, met_name,
             bioMax, bioMean, bioMin,
             gene_assoc, pFBAFlag
            )

  stopifnot(validObject(obj))

  return(obj)
}

#' Methods for Manipulating FBA_greatmod Objects
#'
#' A collection of methods to inspect, modify, and export models defined by the S4 class `FBA_greatmod`.
#' These include functions for adjusting constraints, setting dietary inputs, managing biomass parameters,
#' defining objective functions, extracting exchange reactions, and writing the model to external formats.
#'
#' @seealso \code{\link{FBA_greatmod}} for the class definition, \code{\link{showFBAmethods}} for the list of the methods.
#'
#' @rdname FBA_greatmod-methods
#' @docType methods
#' @return showFBAmethods returns a character vector of method names used for FBA_greatmod objects
#' @aliases showFBAmethods, FBA_greatmod-method
#' @export

showFBAmethods <- function() {
  return(c(
    "setObjFun", "getExchangesR", "getConstraints", "setConstraints",
    "setDiet", "setDiet.name", "writeFBAfile", "setBiomassParameters",
    "setPFbaGeneOption"
  ))
}


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


setMethod(f = "initialize",
          signature("FBA_greatmod"),
          definition = function(.Object, S, ub, lb, obj_fun,
                                react_name=NULL,
                                met_name=NULL,
                                bioMax=-1,
                                bioMean=-1,
                                bioMin=-1,
                                gene_assoc=NULL,
                                pFBAFlag = -1
                                )
          {
            # Campi già esistenti...
            .Object@S       <- as.matrix(S)
            .Object@lowbnd  <- lb
            .Object@uppbnd  <- ub
            .Object@obj_coef<- obj_fun
            .Object@bioMax  <- bioMax
            .Object@bioMean <- bioMean
            .Object@bioMin  <- bioMin
            .Object@pFBAFlag  <- pFBAFlag

            ncol <- ncol(.Object@S)
            nrow <- nrow(.Object@S)

            # Reazioni (react_name)
            if(is.null(react_name)) {
              .Object@react_id = paste0("Reaction", 1:ncol)
            } else {
              .Object@react_id = react_name
            }
            # Metaboliti (met_name)
            if(is.null(met_name)) {
              .Object@met_id = paste0("Metabolite", 1:nrow)
            } else {
              .Object@met_id = met_name
            }

            # NUOVO: assegno lo slot gene_assoc
            if(!is.null(gene_assoc)) {
              if(length(gene_assoc) != ncol) {
                stop("gene_assoc must have length equal to the number of reactions (ncol of S).")
              }
              .Object@gene_assoc <- gene_assoc
            } else {
              # se non lo passi, lo inizializzi a zero o vector numeric(0)
              .Object@gene_assoc <- rep(0, ncol)
            }

            return(.Object)
          }
)


#' Set Objective Function for FBA Model
#'
#' This method assigns the objective coefficients based on the specified reaction names.
#'
#' @aliases setObjFun,FBA_greatmod-method
#' @param theObject An object of class `FBA_greatmod`
#' @param reaction_name A character vector of reaction names to set as objective
#' @return setObjFun returns the updated `FBA_greatmod` object with modified objective coefficients
#' @export
#' @examples
#' \dontrun{
#' model <- setObjFun(model, reaction_name = "EX_glc__D_e")
#' }

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

#' Identify Exchange Reactions
#'
#' getExchangesR returns the list of exchange reactions, i.e., reactions involving only one metabolite,
#' which are typically used for input/output in FBA models.
#'
#' @param theObject An object of class `FBA_greatmod`
#' @return A character vector of reaction IDs classified as exchange reactions.
#' @rdname FBA_greatmod-methods
#' @aliases getExchangesR,FBA_greatmod-method
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

#' Get Constraints for a Specific Reaction
#'
#' Retrieves the lower and upper bounds for a given reaction name in the model.
#'
#' @param theObject An object of class `FBA_greatmod`
#' @param reaction.name A character string identifying the reaction
#' @return A numeric vector of length 2 containing the lower and upper bounds
#' @rdname FBA_greatmod-methods
#' @aliases getConstraints,FBA_greatmod-method
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

#' Set New Constraints for a Specific Reaction
#'
#' Updates the lower and upper bounds for a given reaction in the FBA model.
#'
#' @param theObject An object of class `FBA_greatmod`
#' @param reaction.name A character string identifying the reaction to modify
#' @param newConstraints Numeric vector of length 2 with the new bounds (lower, upper)
#' @return The modified object of class `FBA_greatmod`
#' @rdname FBA_greatmod-methods
#' @aliases setConstraints,FBA_greatmod-method
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

#' Apply Dietary Constraints to FBA Model
#'
#' Updates the reaction bounds in an FBA model based on a provided dietary file or data frame.
#'
#' @param theObject A `FBA_greatmod` object
#' @param dietf Path to a tab-delimited diet file or a data.frame with three columns: reaction ID, lower bound, upper bound.
#' @return A modified `FBA_greatmod` object with updated reaction constraints and recorded diet metadata.
#' @rdname FBA_greatmod-methods
#' @aliases setDiet,FBA_greatmod-method
#' @export
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

#' Manually Set the Diet Name in an FBA Model
#'
#' Updates the internal slot \code{dietName} with a user-defined label.
#'
#' @param theObject A `FBA_greatmod` object
#' @param diet_name A character string representing the diet name
#' @return A modified `FBA_greatmod` object
#' @rdname FBA_greatmod-methods
#' @aliases setDiet.name,FBA_greatmod-method
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


#' Write FBA Model to File
#'
#' Exports the model data to a text file for use with external tools or simulators. Delegates to a compiled C++ function.
#'
#' @param theObject A `FBA_greatmod` object
#' @param fba_fname Filename to save the model (default = "fba_file.txt")
#' @param dest_dir Destination directory for output file
#' @param ... Additional arguments passed to C++ backend
#' @rdname FBA_greatmod-methods
#' @aliases writeFBAfile,FBA_greatmod-method
#' @export
#' @import Rcpp

setGeneric(name="writeFBAfile",
           def=function(theObject,fba_fname,dest_dir,...)
           {
             standardGeneric("writeFBAfile")
           }
)



setMethod(f="writeFBAfile",
  signature=c("FBA_greatmod", "character"),
  definition=function(theObject, fba_fname=NULL, dest_dir=NULL, ...) {
    if (is.null(fba_fname)) {
      fba_fname <- "fba_file.txt"
    }

    # Usa theObject invece di model
    S         <- as.matrix(theObject@S)
    react_id  <- unlist(theObject@react_id)
    obj_coef  <- c(theObject@obj_coef)
    lowbnd    <- as.matrix(theObject@lowbnd)
    uppbnd    <- as.matrix(theObject@uppbnd)
    bioMax    <- theObject@bioMax
    bioMean   <- theObject@bioMean
    bioMin    <- theObject@bioMin
    gen_assoc <- theObject@gene_assoc
    pFBAFlag    <- theObject@pFBAFlag

    # Dichiariamo la variabile model_name (così come serve nel C++):
    model_name <- fba_fname  # se vuoi che il file si chiami come fba_fname
    write <- TRUE

    ncol <- ncol(S)
    nrow <- nrow(S)

    # Creiamo la matrice rb di vincoli
    b <- matrix(0, nrow = ncol, ncol = 1)
    rb <- cbind(b,b)

    Rcpp::sourceCpp(system.file("AuxiliarFunctions/writeFBAfile.cpp", package = "epimodFBAfunctions"))
    # Chiamata alla funzione C++:
    writeModelCpp(
      S          = S,
      react_id   = react_id,
      obj_coef   = obj_coef,
      lowbnd     = lowbnd,
      uppbnd     = uppbnd,
      rb         = rb,
      gene_assoc = gen_assoc,
      model_name = model_name,
      write      = write,
      wd         = dest_dir,
      bioMax     = bioMax,
      bioMean    = bioMean,
      bioMin     = bioMin,
      pFBAFlag   = pFBAFlag
    )
  }
)

#' Set Biomass Summary Values
#'
#' Assigns values for biomass max, mean, and min slots used in model tracking or export.
#'
#' @param object A `FBA_greatmod` object
#' @param bioMax Numeric value for maximum biomass
#' @param bioMean Numeric value for mean biomass
#' @param bioMin Numeric value for minimum biomass
#' @return Modified `FBA_greatmod` object
#' @rdname FBA_greatmod-methods
#' @aliases setBiomassParameters,FBA_greatmod-method
#' @export


setGeneric("setBiomassParameters", function(object, bioMax, bioMean, bioMin) standardGeneric("setBiomassParameters"))

setMethod("setBiomassParameters", signature = "FBA_greatmod", definition = function(object, bioMax, bioMean, bioMin) {
  object@bioMax <- bioMax
  object@bioMean <- bioMean
  object@bioMin <- bioMin
  return(object)
})

#' Set Gene Option for Parsimonious FBA (pFBA)
#'
#' Updates the \code{pFBAFlag} slot in the model, controlling how gene associations are handled in downstream pFBA analysis.
#'
#' @param object A `FBA_greatmod` object
#' @param geneOption A numeric flag for pFBA configuration
#' @return Modified `FBA_greatmod` object
#' @rdname FBA_greatmod-methods
#' @aliases setPFbaGeneOption,FBA_greatmod-method
#' @export

setGeneric("setPFbaGeneOption", function(object, geneOption) standardGeneric("setPFbaGeneOption"))

setMethod("setPFbaGeneOption", signature = "FBA_greatmod", definition = function(object, geneOption) {
  object@pFBAFlag <- geneOption
  return(object)
})
