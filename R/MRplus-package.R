#' @keywords internal
#' @importFrom dplyr ungroup %>% mutate summarise group_by
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_text
#'   scale_size_continuous scale_color_gradientn scale_y_log10
#'   labs theme_classic theme element_text element_rect element_line element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
#' @importFrom scales comma
#' @importFrom TwoSampleMR extract_instruments extract_outcome_data harmonise_data mr
"_PACKAGE"
