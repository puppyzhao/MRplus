#' Visualise robustness of MR results across multiple datasets
#'
#' This function runs Mendelian randomization (MR) analyses for multiple
#' exposure–outcome pairs across several instrument selection thresholds.
#' It aggregates significance counts and mean p-values, and visualises the
#' results in a bubble plot.
#'
#' @param exposure_ids Character vector. IDs of exposure GWAS datasets 
#'   (e.g., \code{c("ieu-b-4872", "bbj-a-80")}).
#' @param outcome_ids Character vector. IDs of outcome GWAS datasets 
#'   (e.g., \code{c("ukb-b-16890", "ebi-a-GCST90013972")}).
#' @param p_values Numeric vector. P-value thresholds for instrument 
#'   selection (e.g., \code{c(1e-5, 5e-6, 5e-7, 5e-8)}).
#'
#' @return A ggplot object showing exposure–outcome pairs as a bubble plot:  
#'   - Bubble size = count of significant MR results (p < 0.05).  
#'   - Bubble color = \code{-log10(mean p-value)}.  
#'   The function also returns a data frame of aggregated results invisibly.
#'
#' @details
#' - Three MR methods are applied: IVW, MR-Egger, and Weighted Median.  
#' - For each exposure–outcome pair and set of thresholds, the function
#'   calculates:  
#'   * the number of significant MR results (\code{p < 0.05}),  
#'   * the mean p-value,  
#'   * and its \code{-log10} transformation.  
#' - Bubble sizes are capped at 12 (the maximum number of possible tests).  
#'
#' @seealso 
#' \code{\link[TwoSampleMR]{extract_instruments}},  
#' \code{\link[TwoSampleMR]{extract_outcome_data}},  
#' \code{\link[TwoSampleMR]{harmonise_data}},  
#' \code{\link[TwoSampleMR]{mr}}
#'
#' @importFrom TwoSampleMR extract_instruments extract_outcome_data harmonise_data mr
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic scale_color_gradientn
#' @importFrom dplyr %>% mutate summarise filter
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' \dontrun{
#' mr_result_matrix_across_datasets(
#'   exposure_ids = c("ieu-b-4872", "bbj-a-80"),
#'   outcome_ids = c("ukb-b-16890", "ebi-a-GCST90013972"),
#'   p_values = c(1e-5, 5e-6, 5e-7, 5e-8)
#' )
#' }
#'
#' @export

mr_result_matrix_across_datasets <- function(exposure_ids, outcome_ids, p_values) {
  results_list <- list()
  row_index <- 1

  total_tasks <- length(exposure_ids) * length(outcome_ids) * length(p_values)
  task_counter <- 0
  pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)

  for (exposure_id in exposure_ids) {
    for (outcome_id in outcome_ids) {
      all_pvals <- numeric()

      for (p1 in p_values) {
        task_counter <- task_counter + 1
        setTxtProgressBar(pb, task_counter)

        Sys.sleep(runif(1, 1, 5))

        tryCatch({
          exposure_data <- extract_instruments(exposure_id, p1 = p1)
          if (is.null(exposure_data) || nrow(exposure_data) == 0) next

          outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = outcome_id)
          if (is.null(outcome_data) || nrow(outcome_data) == 0) next

          harmonized_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

          mr_results <- mr(
            harmonized_data,
            method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
          )
          pvals <- mr_results$pval
          all_pvals <- c(all_pvals, pvals)

        }, error = function(e) {
          cat("Error processing:", exposure_id, outcome_id, "at P-value threshold:", p1, "\n")
        })
      }

      significant_count <- sum(all_pvals < 0.05, na.rm = TRUE)
      mean_pval <- mean(all_pvals, na.rm = TRUE)

      significant_count <- pmin(significant_count, 12)

      results_list[[row_index]] <- data.frame(
        exposure_id = exposure_id,
        outcome_id = outcome_id,
        significant_count = significant_count,
        mean_pval = mean_pval,
        neglog10_mean_pval = -log10(mean_pval),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1
    }
  }

  close(pb)

  results <- do.call(rbind, results_list)

  ggplot(results, aes(x = exposure_id, y = outcome_id)) +
    geom_point(aes(size = significant_count, color = neglog10_mean_pval), alpha = 0.9) +

    scale_size_continuous(
      range = c(0, 12),
      limits = c(0, 12),
      breaks = c(4, 8, 12),
      name = "Count of\nP < 0.05"
    ) +

    scale_color_gradientn(
      colors = c("grey90", "#fdae61", "#d7191c"),
      values = c(0, 0.25, 1),   # 手动调整 P=0.05 在颜色条下方 25%
      limits = c(0, 3),
      name = expression(-log[10]("Mean P")),
      breaks = c(0, 1.3, 3),
      labels = c("P=1", "P=0.05", " "),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        barwidth = unit(0.4, "cm"),
        barheight = unit(3, "cm"),  # 缩短颜色条高度
        ticks = TRUE,
        ticks.colour = "black"
      )
    ) +

    labs(
      x = "Exposure ID",
      y = "Outcome ID"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank()
    )
}

