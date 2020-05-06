#' Relevel the input data frame; useful for post-filtering organization
#' 
#' @param x The input data frame
#' 
#' @return Releveled data frame
relevel_df <- function(x) {
    for (cn in colnames(x)) {
        if (class(x[[cn]]) == "factor") {
            x[[cn]] <- as.factor(as.character(x[[cn]]))
        }
    }
    return(x)
}

#' Transform input data given covariates and surrogate variables
#' 
#' @param cts Input data
#' @param mm Input model matrix
#' @param sv Input surrogate variable matrix
#'
#' @return Transformed counts matrix
transform_data <- function(cts, mm, sv) {
    # Operating under the assumption that mm is completely organized
    sample_names <- rownames(mm)
    cts <- cts[, sample_names]
    X <- cbind(mm, sv[sample_names, ])
    Beta <- (solve(t(X) %*% X) %*% t(X)) %*% t(cts)
    return(cts - t(as.matrix(X[, -c(1:ncol(mm))]) %*% Beta[-c(1:ncol(mm)), ]))
}