#' Compare MR results across different IV selection thresholds
#'
#' This function performs Mendelian randomization (MR) analyses for a given
#' exposure–outcome pair across multiple instrument selection thresholds.
#' It extracts instruments at each p-value threshold, harmonises exposure–
#' outcome data, runs MR with different methods (IVW, MR-Egger, and Weighted
#' Median), and visualises the resulting p-values.
#'
#' @param exposure_id Character. The ID of the exposure GWAS 
#'   (e.g., "ieu-b-4872").
#' @param outcome_id Character. The ID of the outcome GWAS 
#'   (e.g., "ukb-b-16890").
#'
#' @return A ggplot object showing MR p-values across different IV selection
#'   thresholds and methods.  
#'   The function also prints the underlying results table to the console.
#'
#' @details 
#' - Tested thresholds are fixed at \code{c(1e-5, 5e-6, 5e-7, 5e-8)}.  
#' - Three MR methods are included:  
#'   * IVW (\code{mr_ivw})  
#'   * MR-Egger regression (\code{mr_egger_regression})  
#'   * Weighted Median (\code{mr_weighted_median})  
#' - P-values are plotted on a \code{log10} scale (minimum set to 1e-4).  
#' - Combinations that fail are imputed with \code{p = 1}.  
#'
#' @seealso 
#' \code{\link[TwoSampleMR]{extract_instruments}},  
#' \code{\link[TwoSampleMR]{extract_outcome_data}},  
#' \code{\link[TwoSampleMR]{harmonise_data}},  
#' \code{\link[TwoSampleMR]{mr}}
#'
#' @importFrom TwoSampleMR extract_instruments extract_outcome_data harmonise_data mr
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_classic
#' @importFrom dplyr %>% mutate filter
#'
#' @examples
#' \dontrun{
#' mr_result_across_pval_thresholds("ieu-b-4872", "ukb-b-16890")
#' }
#'
#' @export


mr_result_across_pval_thresholds <- function(exposure_id, outcome_id) {

  # Define outcome(s) and IV thresholds
  outcome_list <- c(outcome_id)
  p_values <- c(1e-5, 5e-6, 5e-7, 5e-8)

  # Initialize result data frame
  ivw_pvals <- data.frame(
    p_threshold = numeric(),
    method = character(),
    pval = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over p-value thresholds
  for (p1 in p_values) {
    cat("\nRunning with P-value threshold:", p1, "\n")

    exposure_list <- c(exposure_id)
    for (exposure0 in exposure_list) {
      tryCatch({
        # Step 1: Extract instruments
        exposure6 <- extract_instruments(exposure0, p1 = p1)

        if (is.null(exposure6) || nrow(exposure6) == 0) {
          cat("No SNPs found for Exposure:", exposure0, "at P-value threshold:", p1, "\n")
          next
        }

        # Step 2: Loop through outcomes
        for (outcome0 in outcome_list) {
          tryCatch({
            # Step 3: Extract outcome data
            outcome1 <- extract_outcome_data(snps = exposure6$SNP, outcomes = outcome0)
            if (is.null(outcome1) || nrow(outcome1) == 0) {
              cat("No outcome data found for Exposure:", exposure0, "Outcome:", outcome0, "\n")
              next
            }

            # Step 4: Harmonize data
            dat2 <- harmonise_data(exposure_dat = exposure6, outcome_dat = outcome1)

            # Step 5: Run MR using different methods
            mr_result_ivw <- mr(dat2, method_list = c("mr_ivw"))
            mr_result_egger <- mr(dat2, method_list = c("mr_egger_regression"))
            mr_result_weighted_median <- mr(dat2, method_list = c("mr_weighted_median"))

            # Step 6: Extract p-values
            ivw_pval <- mr_result_ivw$pval[1]
            egger_pval <- mr_result_egger$pval[1]
            weighted_median_pval <- mr_result_weighted_median$pval[1]

            # Step 7: Save results
            ivw_pvals <- rbind(ivw_pvals, data.frame(p_threshold = p1, method = "IVW", pval = ivw_pval))
            ivw_pvals <- rbind(ivw_pvals, data.frame(p_threshold = p1, method = "Egger", pval = egger_pval))
            ivw_pvals <- rbind(ivw_pvals, data.frame(p_threshold = p1, method = "Weighted Median", pval = weighted_median_pval))

          }, error = function(e) {
            cat("Error processing Exposure:", exposure0, "Outcome:", outcome0, "\n")
            cat("Error message:", e$message, "\n")
          })
        }
      }, error = function(e) {
        cat("Error processing Exposure:", exposure0, "at P-value threshold:", p1, "\n")
        cat("Error message:", e$message, "\n")
      })
    }
  }

  # Fill missing combinations with p = 1
  for (p1 in p_values) {
    for (method in c("IVW", "Egger", "Weighted Median")) {
      if (!any(ivw_pvals$p_threshold == p1 & ivw_pvals$method == method)) {
        ivw_pvals <- rbind(ivw_pvals, data.frame(p_threshold = p1, method = method, pval = 1))
      }
    }
  }

  # Format and clean
  ivw_pvals$p_threshold <- factor(ivw_pvals$p_threshold, levels = rev(p_values))
  ivw_pvals$pval[ivw_pvals$pval < 0.0001] <- 0.0001
  ivw_pvals$pval[is.na(ivw_pvals$pval)] <- 1

  print(ivw_pvals)

  # Plot: MR P-values across different thresholds and methods
  ggplot(ivw_pvals, aes(x = p_threshold, y = pval, color = method, shape = method, group = method)) +
    geom_line(size = 1.2, alpha = 0.7) +
    geom_point(size = 4, alpha = 0.7) +
    scale_y_log10(
      limits = c(0.0001, 1),
      breaks = c(0.001, 0.01, 0.05, 1),
      labels = label_number(scale = 1, accuracy = 0.001)
    ) +
    scale_color_manual(
      values = c(
        "IVW" = "#264653",
        "Egger" = "#6a4c93",
        "Weighted Median" = "#ffba08"
      )
    ) +
    scale_shape_manual(
      values = c(
        "IVW" = 16,
        "Egger" = 17,
        "Weighted Median" = 15
      )
    ) +
    labs(
      x = "IV Selection Threshold",
      y = "MR P-value",
      color = "Method",
      shape = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x = element_text(size = 14, family = "Helvetica"),
      axis.text.y = element_text(size = 14, family = "Helvetica"),
      axis.title = element_text(size = 14, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 14, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 12, family = "Helvetica"),
      legend.position = "right",
      legend.box.spacing = unit(0.3, "cm"),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_line(color = "gray95", size = 0.25),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    geom_hline(yintercept = 0.05, color = "gray50", linetype = "dashed", size = 1)
}
