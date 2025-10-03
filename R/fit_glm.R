#' Title: Fitting generalized linear models for the best model
#'
#' Description: `glm.best` is used to fit generalized linear model for the best model provided by `modelselect.glm`.
#'
#' @param object the model selection result from `modelselect.glm`.
#' @param family a character string naming a family function describing the error distribution to be used in the model.
#' @param method the criteria to do model select.
#'               `method = "models"` selects the best model by the highest posterior probabilities.
#'               `method = "variables"` selects the variables in the best model by the posterior inclusion probabilities which are larger than the threshold.
#' @param threshold The threshold for variable selection. The variables with posterior inclusion probability larger than the threshold are selected in the best model. The default is 0.95.
#' @param x,y logicals. If `TRUE` the corresponding components (the best model predictor matrix, the response) of the fit are returned.
#' @return An object of class `"glm"`, which is a list containing the following components:
#' \describe{
#'   \item{\code{coefficients}}{a named vector of coefficients.}
#'   \item{\code{residuals}}{the working residuals, that is the residuals in the final iteration of the IWLS fit.}
#'   \item{\code{fitted.values}}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'   \item{\code{rank}}{the numeric rank of the fitted linear model.}
#'   \item{\code{family}}{the family object used.}
#'   \item{\code{linear.predictors}}{the linear fit on the link scale.}
#'   \item{\code{deviance}}{up to a constant, minus twice the maximized log-likelihood.}
#'   \item{\code{aic}}{A version of Akaike's An Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters, computed by the \code{aic} component of the family.}
#'   \item{\code{null.deviance}}{The deviance for the null model, comparable with \code{deviance}. The null model will include the offset, and an intercept if there is one in the model.}
#'   \item{\code{iter}}{the number of iterations of IWLS used.}
#'   \item{\code{weights}}{the working weights, that is the weights in the final iteration of the IWLS fit.}
#'   \item{\code{prior.weights}}{the weights initially supplied, a vector of 1s if none were.}
#'   \item{\code{df.residual}}{the residual degrees of freedom.}
#'   \item{\code{df.null}}{the residual degrees of freedom for the null model.}
#'   \item{\code{y}}{if requested, the response vector used.}
#'   \item{\code{converged}}{logical. Was the IWLS algorithm judged to have converged?}
#'   \item{\code{boundary}}{logical. Is the fitted value on the boundary of the allowable values?}
#'   \item{\code{model}}{if requested (the default), the model frame used.}
#'   \item{\code{call}}{the matched call.}
#'   \item{\code{formula}}{the formula supplied.}
#'   \item{\code{terms}}{the \code{\link{terms.object}} used.}
#'   \item{\code{data}}{the data argument.}
#'   \item{\code{threshold}}{the threshold used for method = "variables".}
#' }

#' @export
glm.best <- function(object, family, method = "models", threshold = 0.95, x = FALSE, y = FALSE){

  if(is.null(object)){
    stop("please provide object.")
  }

  if(is.null(family)){
    stop("please provide family.")
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
    fit <- stats::glm(formula = stats::as.formula(paste0(response_name,"~1")),
                      data = data, family=family, x=x, y=y)
    fit$call <- call("glm", formula = stats::as.formula(paste0(response_name, "~1")), family = family)
  }else{
    fit <- stats::glm(formula = stats::as.formula(paste0(response_name,"~", paste(fixedeffect_name, collapse ="+"))),
                      data = data, family=family, x=x, y=y)
    fit$call <- call("glm", formula = stats::as.formula(paste0(response_name,"~", paste(fixedeffect_name, collapse ="+"))), family = family)
  }

  cat("Call:\n")
  print(fit$call)

  cat("\nCoefficients:\n")
  print(fit$coefficients)
  invisible(fit)

}


