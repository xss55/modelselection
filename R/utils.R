#' Function for checking model hierarchy
#' @return model hierarchy
#' @keywords internal


check_model_hierarchy <- function(model) {
    terms_obj <- stats::terms(model)
    term_labels <- attr(terms_obj, "term.labels")

    interactions <- term_labels[grepl(":", term_labels)]

    issues_found <- FALSE

    for (interaction in interactions) {
      parts <- unlist(strsplit(interaction, ":"))
      n_parts <- length(parts)

      for (k in 1:(n_parts - 1)) {
        combos <- utils::combn(parts, k, simplify = FALSE)
        for (combo in combos) {
          lower_term <- paste(combo, collapse = ":")
          if (!(lower_term %in% term_labels)) {
            issues_found <- TRUE
          }
        }
      }
    }
    return(issues_found)
}
