#' Convert expression matrix to z-scores
#'
#' @param expr Expression matrix
#' 
#' @return Converted expression matrix
expr_to_z <- function(expr) {
	require(magrittr)
	expr %>%
		sweep(., 2, colMeans(expr), "-") %>%
		sweep(., 2, apply(expr, 2, sd), "/")
}

#' Return developmental period divisions for plotting
#'
#' @param metadata Metadata for period definitions 
#'                 (should include "Days" and "Period" fields)
#' 
#' @return Numeric vector of period breakpoints
get_period_breaks <- function(metadata) {
	sapply(
    	2:length(unique(metadata$Period)),
    	function(i) {
	        p1 <- sort(unique(metadata$Period))[[i - 1]]
	        p2 <- sort(unique(metadata$Period))[[i]]
	        max_p1 <- max(pull(filter(metadata, Period == p1), Days))
	        min_p2 <- min(pull(filter(metadata, Period == p2), Days))
	        mean(c(max_p1, min_p2))
	    }
	)
}

#' Add standard geoms/parts to an expression profile plot
#'
#' @param plt Input plot
#' @param metadata Metadata for period definitions 
#'                 (should include "Days" and "Period" fields)
#' 
#' @return Filled expression profile ggplot object
expression_plot_fill <- function(plt, metadata) {
	plt +
		geom_vline(
			xintercept = get_period_breaks(metadata), colour = "black"
		) +
		annotate(
			geom = "rect",
			xmin = min(metadata$Days), 
			xmax = get_period_breaks(metadata)[[6]], 
			ymin = -Inf, ymax = Inf, 
			fill = "red", alpha = 0.1
		) +
		annotate(
			geom = "rect",
			xmin = get_period_breaks(metadata)[[6]], 
			xmax = max(metadata$Days), 
			ymin = -Inf, ymax = Inf, 
			fill = "blue", alpha = 0.1
		) +
		scale_x_continuous(
			trans = "log1p",
			breaks = c(50, 100, 1000, 10000),
			expand = c(0, 0)
		)
}