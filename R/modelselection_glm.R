#' Title: Model selection for generalized linear models
#'
#' Description: use BIC to do model selection.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#'                A typical model has the form `response ~ terms` where response is the (numeric) `response` vector and terms is a series of terms which specifies a linear predictor for `response`.
#'                A terms specification of the form `first + second` indicates all the terms in `first` together with all the terms in `second` with duplicates removed.
#'                A specification of the form `first:second` indicates the set of terms obtained by taking the interactions of all terms in `first` with all terms in `second`.
#'                The specification `first*second` indicates the cross of `first` and `second.` This is the same as `first + second + first:second`.
#' @param data  an data frame containing the variables in the model.
#' @param family  a character string naming a family function describing the error distribution to be used in the model.
#' @param GA_var if the number of variables is smaller than `GA_var`, then do exhaustive model search, otherwise use genetic algorithm to do stochastic model search.
#' @param maxiterations  the maximum number of iterations to run before the GA search is halted.
#' @param runs_til_stop 	the number of consecutive generations without any improvement in the best fitness value before the GA is stopped.
#' @param monitor 	a logical defaulting to TRUE showing the evolution of the search. If monitor = FALSE, any output is suppressed.
#' @param popSize the population size.
#' @return `modelselect.glm` returns a list containing the following components:
#'
#' \describe{
#'   \item{\code{models}}{A data frame of candidate models' BIC and posterior probabilities, sorted by decreasing posterior probability}
#'   \item{\code{variables}}{A data frame of candidate variables' posterior inclusion probabilities}
#'   \item{\code{data}}{The data with variables in the formula.}
#' }
#'
#' The function `glm.best` is used to obtain the linear fitting to the best model by posterior probability or by controlling variables' posterior inclusion probabilities.
#'
#' @export

modelselect.glm <- function(formula, data, family, GA_var=16, maxiterations = 2000, runs_til_stop = 1000, monitor = TRUE, popSize = 100){

  if(is.null(formula)){
    stop("Please provide formula.")
  }

  if(is.null(data)){
    stop("Please provide data.")
  }

  data <- as.data.frame(data)
  terms_obj <- stats::terms(formula, data = data)
  terms_in_formula <- all.vars(terms_obj)
  missing_vars <- setdiff(terms_in_formula, names(data))

  if (length(missing_vars) > 0) {
    stop(paste("The following variables are missing from the data:", paste(missing_vars, collapse = ", ")))
  }

  data <- data[, terms_in_formula]
  if(sum(is.na(data))){
    stop("There is missing data.")
  }

  response_name <- as.character(attr(terms_obj, "variables"))[attr(terms_obj, "response") + 1]

  X_colname <- attr(terms_obj, "term.labels")

  K <- length(X_colname)
  if(K < GA_var){
    binary_matrix <- rep(list(0:1), K)
    binary_matrix <- as.matrix(expand.grid(binary_matrix))

    model_bic_prob <- c()
    for (i in 0:(2^K-1)){
      if(i == 0){
        fit <- stats::glm(formula = paste0(response_name, "~1"), data = data, family = family)
        bic <-  stats::BIC(fit)
        if(is.na(bic)){
          model_bic_prob <- rbind(model_bic_prob, c(binary_matrix[i+1,], Inf))
        }else{
          model_bic_prob <- rbind(model_bic_prob, c(binary_matrix[i+1,], bic))
        }

      }else{

        sub_colname <- X_colname[which(binary_matrix[i+1,]==1)]
        fit <- stats::glm(formula = paste0(paste0(response_name, "~"), paste(sub_colname,collapse = "+")), data = data, family = family)

        if(!check_model_hierarchy(fit)){
          bic <-  stats::BIC(fit)
          if(is.na(bic)){
            model_bic_prob <- rbind(model_bic_prob, c(binary_matrix[i+1,], Inf))
          }else{
            model_bic_prob <- rbind(model_bic_prob, c(binary_matrix[i+1,], bic))
          }

        }
      }
    }

    bic <- model_bic_prob[,ncol(model_bic_prob)]
    bestmodel <- model_bic_prob[which.min(bic),-ncol(model_bic_prob)]

    dat <- model_bic_prob[,1:K]

  }else{

    fitness_ftn <- function(string){
      if(sum(string) == 0){
        return(-stats::BIC(stats::glm(formula = paste0(response_name, "~1"), data = data, family = family)))
      }else{

        model <- which(string==1)
        k <- length(model)

        sub_colname <- X_colname[model]
        fit <-stats::glm(formula = paste0(paste0(response_name, "~"), paste(sub_colname,collapse = "+")), data = data, family = family)
        if(!check_model_hierarchy(fit)){

          bic <-  stats::BIC(fit)
          if(is.na(bic)){
            return(-Inf)
          }else{
            return(-bic)
          }

        }
      }
    }


    if(K > 99){
      suggestedsol <- diag(K)
      tmp_BIC <- vector()
      for(i in 1:K){
        model <- which(suggestedsol[i,]==1)
        k <- length(model)

        sub_colname <- X_colname[model]
        tmp_BIC[i] <- -stats::BIC(stats::glm(formula = paste0(paste0(response_name, "~"), paste(sub_colname,collapse = "+")), data = data, family = family))
      }
      suggestedsol <- rbind(0,suggestedsol[order(tmp_BIC,decreasing = TRUE)[1:99],])
    }else{
      suggestedsol <- rbind(0,diag(K))
    }


    fitness_ftn <- memoise::memoise(fitness_ftn)
    ans <- GA::ga("binary", fitness = fitness_ftn,
                  nBits = K,maxiter = maxiterations,popSize = popSize,
                  elitism = min(c(10,2^K)),run = runs_til_stop,suggestions = suggestedsol,monitor = monitor)
    memoise::forget(fitness_ftn)
    dat <- ans@population
    dupes <- duplicated(dat)
    dat <- dat[!dupes,]
    ans@fitness <- ans@fitness[!dupes]
    bic <- -ans@fitness
    model_bic_prob <- cbind(dat,bic)
    model_bic_prob <- stats::na.omit(model_bic_prob)
    bic <- model_bic_prob[,ncol(model_bic_prob)]
    dat <- model_bic_prob[,1:(ncol(model_bic_prob)-1)]
    bestmodel <- dat[which.min(bic),]

  }

  postprob_temp <- exp(-0.5*(bic-min(bic)))
  postprob <- postprob_temp/sum(postprob_temp)
  model_bic_prob <- cbind(model_bic_prob, postprob)
  model_bic_prob <- model_bic_prob[!(bic == Inf),]

  var_pip <- X_colname
  if(K == 1){
    pip <- sum(dat * postprob)
  }else{
    pip <- colSums(dat * postprob)
  }
  var_pip <- data.frame("Var"=var_pip, "PIP" = round(pip,3))

  model_bic_prob <- as.data.frame(model_bic_prob)
  if(is.null(X_colname)){
    colnames(model_bic_prob) <- c(paste("X",1:K,sep=""),"BIC","MPP")
  }else{
    colnames(model_bic_prob) <- c(X_colname,"BIC","MPP")

  }
  model_bic_prob <- model_bic_prob[order(model_bic_prob$MPP, decreasing = T),]
  model_bic_prob[,"MPP"] <- round(model_bic_prob[,"MPP"], 3)
  rownames(model_bic_prob) <- 1:nrow(model_bic_prob)

  attr(model_bic_prob, "id") <- "models"
  attr(var_pip, "id") <- "variables"
  attr(data, "id") <- response_name

  cat("The Best Model:\n", colnames(model_bic_prob)[1:K][model_bic_prob[1,1:K] == 1], "\n")
  cat("BIC:\n", model_bic_prob[1,K+1], "\n")
  cat("MPP:\n", model_bic_prob[1,K+2])

  invisible(list("models" = model_bic_prob, "variables" = var_pip, "data" = data))

}



