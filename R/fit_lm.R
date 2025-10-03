#' Title: Fitting linear models for the best model
#'
#' Description: `lm.best` is used to fit linear model for the best model provided by `modelselect.lm`.
#'
#' @param object the model selection result from `modelselect.lm`.
#' @param method the criteria to do model select.
#'               `method = "models"` selects the best model by the highest posterior probabilities.
#'               `method = "variables"` selects the variables in the best model by the posterior inclusion probabilities which are larger than the threshold.
#' @param threshold The threshold for variable selection. The variables with posterior inclusion probability larger than the threshold are selected in the best model. The default is 0.95.
#' @param x,y logicals. If `TRUE` the corresponding components (the best model predictor matrix, the response) of the fit are returned.
#' @return An object of class `"lm"`, which is a list containing the following components:
#' \describe{
#'   \item{\code{coefficients}}{A named vector of coefficients.}
#'   \item{\code{residuals}}{The residuals, that is the response minus the fitted values.}
#'   \item{\code{fitted.values}}{The fitted mean values.}
#'   \item{\code{rank}}{The numeric rank of the fitted linear model.}
#'   \item{\code{df.residual}}{The residual degrees of freedom.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{terms}}{The `terms` object used.}
#'   \item{\code{model}}{(If requested) the model frame used.}
#'   \item{\code{qr}}{(If requested) the QR decomposition of the design matrix.}
#'   \item{\code{xlevels}}{(If the model formula includes factors) a record of the levels of the factors.}
#'   \item{\code{contrasts}}{(If the model formula includes factors) the contrasts used.}
#'   \item{\code{offset}}{The offset used.}
#'   \item{\code{threshold}}{the threshold used for method = "variables".}
#' }

#' @export
lm.best <- function(object, method = "models", threshold = 0.95, x = FALSE, y = FALSE){

  if(is.null(object)){
    stop("please provide object.")
  }
  data <- object$data
  n <- nrow(data)
  response_name <- attr(data, "id")

  if(method == "models"){
    model_bic_prob <- object[[method]]
    fixedeffect_name <- colnames(model_bic_prob)
    fixedeffect_name <- fixedeffect_name[which(model_bic_prob[1,1:(length(fixedeffect_name)-2)] ==1)]
  }

  if(method == "variables"){
    var_pip <- object[[method]]
    fixedeffect_name <- var_pip$Var[var_pip$PIP>threshold]
  }

  if(length(fixedeffect_name) == 0){
    fit <- stats::lm(formula = stats::as.formula(paste0(response_name,"~1")),
                     data = data, x=x, y=y)
    fit$call <- call("lm", stats::as.formula(paste0(response_name, "~1")))
  }else{
    fit <- stats::lm(formula = stats::as.formula(paste0(response_name,"~", paste(fixedeffect_name, collapse ="+"))),
                     data = data, x=x, y=y)
    fit$call <- call("lm", stats::as.formula(paste0(response_name,"~", paste(fixedeffect_name, collapse ="+"))))
  }


  cat("Call:\n")
  print(fit$call)

  cat("\nCoefficients:\n")
  print(fit$coefficients)
  invisible(fit)

}


