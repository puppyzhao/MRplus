#' Visualize the number of instruments across LD pruning and p-value thresholds
#'
#' This function extracts instruments for a given exposure from MR-Base
#' (via \code{TwoSampleMR::extract_instruments}) under multiple 
#' p-value thresholds and LD pruning settings. It summarizes the SNP counts 
#' and generates a visualization of how instrument selection is influenced 
#' by these parameters.
#'
#' @param exposure_id Character. The ID of the exposure GWAS (e.g., "ieu-b-4872").
#' @param p_values Numeric vector. A set of p-value thresholds to test for instrument selection.
#'   Default is \code{c(5e-8, 5e-7, 5e-6, 1e-5)}.
#' @param r2_kb_combinations List of numeric vectors. Each element should be of the form 
#'   \code{c(r2, kb)} specifying the LD pruning threshold (r2) and window size (kb).
#'   Default is \code{list(c(0.0001, 1000), c(0.001, 10000), c(0.01, 100000))}.
#'
#' @return A ggplot object showing SNP counts under different thresholds and LD settings.  
#'   The function also prints progress information to the console.
#'
#' @details 
#' - The "Default" setting corresponds to \code{r2 = 0.001, kb = 10000}, which is 
#'   the default in \code{TwoSampleMR}.
#' - The output plot shows:  
#'   * Line + points: SNP counts under the default LD pruning  
#'   * Triangles and squares: SNP counts under min and max LD pruning  
#'   * Red error bars: range (min–max) across LD pruning at each p-value threshold  
#'   * Labels: annotated min and max values  
#'
#' @seealso \code{\link[TwoSampleMR]{extract_instruments}}
#'
#' @importFrom TwoSampleMR extract_instruments
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @importFrom dplyr %>% group_by summarise arrange
#'
#' @examples
#' \dontrun{
#' mr_iv_count_across_ld_and_pval_threshold("ieu-b-4872")
#' }
#'
#' @export

mr_iv_count_across_ld_and_pval_threshold <- function(
    exposure_id,
    p_values = c(5e-8, 5e-7, 5e-6, 1e-5),
    r2_kb_combinations = list(
      c(0.0001, 1000),
      c(0.001, 10000),
      c(0.01, 100000)
    )
) {
  # Initialize result dataframe
  results <- data.frame(
    p_value = numeric(),
    r2_kb = character(),
    snp_count = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over p-value thresholds and LD pruning parameters
  for (p1 in p_values) {
    cat("Processing p-value threshold:", p1, "\n")
    for (r2_kb in r2_kb_combinations) {
      r2 <- r2_kb[1]
      kb <- r2_kb[2]
      cat("Processing r2 =", r2, ", kb =", kb, "\n")

      tryCatch({
        dat <- extract_instruments(exposure_id, r2 = r2, kb = kb, p1 = p1)
        n <- ifelse(is.null(dat), 0, nrow(dat))
        results <- rbind(
          results,
          data.frame(p_value = p1, r2_kb = paste("r2 =", r2, "kb =", kb), snp_count = n)
        )
      }, error = function(e) {
        cat("Error for p =", p1, ", r2 =", r2, ", kb =", kb, "\n")
        results <- rbind(
          results,
          data.frame(p_value = p1, r2_kb = paste("r2 =", r2, "kb =", kb), snp_count = 0)
        )
      })
    }
  }

  # Summarize results by p-value threshold
  results$p_value <- as.numeric(as.character(results$p_value))

  summary_data <- results %>%
    group_by(p_value) %>%
    summarise(
      default_snp_count = mean(snp_count),
      min_snp_count = min(snp_count),
      max_snp_count = max(snp_count)
    ) %>%
    ungroup() %>%
    mutate(
      p_value = factor(p_value, levels = p_values),
      x_position = as.numeric(p_value)
    )

  # Convert to long format for plotting
  summary_data_long <- summary_data %>%
    pivot_longer(
      cols = c(default_snp_count, min_snp_count, max_snp_count),
      names_to = "type",
      values_to = "snp_count"
    ) %>%
    mutate(
      type = factor(
        case_when(
          type == "default_snp_count" ~ "Default (r2=0.001, kb=10000)",
          type == "min_snp_count" ~ "Min (r2=0.0001, kb=1000)",
          type == "max_snp_count" ~ "Max (r2=0.01, kb=100000)"
        ),
        levels = c(
          "Min (r2=0.0001, kb=1000)",
          "Default (r2=0.001, kb=10000)",
          "Max (r2=0.01, kb=100000)"
        )
      )
    )

  # Plot: Number of SNPs selected under different thresholds and LD pruning settings
  ggplot() +
    geom_line(
      data = subset(summary_data_long, type == "Default (r2=0.001, kb=10000)"),
      aes(x = x_position, y = snp_count, color = type),
      linetype = "solid", size = 1.1, alpha = 0.8
    ) +
    geom_point(
      data = subset(summary_data_long, type == "Default (r2=0.001, kb=10000)"),
      aes(x = x_position, y = snp_count, color = type),
      size = 4, shape = 19, alpha = 0.8
    ) +
    geom_point(
      data = subset(summary_data_long, type == "Min (r2=0.0001, kb=1000)"),
      aes(x = x_position, y = snp_count, color = type),
      size = 4, shape = 17, alpha = 0.7
    ) +
    geom_point(
      data = subset(summary_data_long, type == "Max (r2=0.01, kb=100000)"),
      aes(x = x_position, y = snp_count, color = type),
      size = 4, shape = 15, alpha = 0.7
    ) +
    geom_errorbar(
      data = summary_data,
      aes(x = x_position, ymin = min_snp_count, ymax = max_snp_count),
      width = 0.15, color = "#D32F2F", size = 0.8, alpha = 0.8
    ) +
    geom_text_repel(
      data = summary_data,
      aes(x = x_position + 0.1, y = max_snp_count, label = sprintf("Max: %.0f", max_snp_count)),
      vjust = -1, color = "#D32F2F", size = 4
    ) +
    geom_text_repel(
      data = summary_data,
      aes(x = x_position - 0.1, y = min_snp_count, label = sprintf("Min: %.0f", min_snp_count)),
      vjust = 2, color = "#D32F2F", size = 4
    ) +
    scale_x_continuous(
      name = "IV selection threshold",
      breaks = summary_data$x_position,
      labels = format(p_values, scientific = TRUE),
      limits = c(0.5, length(p_values) + 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    ylab("Number of selected SNPs") +
    scale_color_manual(
      name = "LD pruning setting",
      values = c(
        "Min (r2=0.0001, kb=1000)" = "#D32F2F",
        "Default (r2=0.001, kb=10000)" = "#005A88",
        "Max (r2=0.01, kb=100000)" = "#A00000"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(size = 0.4, linetype = "dashed", color = "grey80"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    coord_cartesian(clip = "off")
}
