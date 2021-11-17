#' plot_peak_saturation_lines
#'
#' @param cnt_dt data.table returned by load_counts()
#' @param color_by variable in cnt_dt to map to color.
#'
#' @return ggplot of peak count vs read count.
#' @export
#' @import ggplot2
#' @examples
plot_peak_saturation_lines = function(cnt_dt, color_by = NULL){
  cnt_dt$group = paste(cnt_dt$out_dir, cnt_dt$sample, cnt_dt$peak_stat)
  lab_dt = data.table::rbindlist(lapply(split(cnt_dt, cnt_dt$group), function(x){
    x[x$read_count == max(x$read_count),]
  }))
  lab_dt$read_count_lab = lab_dt$read_count + 8e5
  if(is.null(color_by)){
    p_lines = ggplot(cnt_dt, aes(x = read_count, y = peak_count, group = group)) +
      geom_path() +
      geom_point(data = lab_dt) +
      geom_text(data = lab_dt, aes(label = sample, x = read_count_lab), hjust = 0) +
      scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
      scale_y_continuous(labels = function(x)x/1e3) +
      labs(x = "reads (M)", y = "peaks (k)") +
      cowplot::theme_cowplot() +
      expand_limits(y = 0)
  }else{
    if(!color_by %in% colnames(cnt_dt)) stop(color_by, " not in colnames of cnt_dt.")
    p_lines = ggplot(cnt_dt, aes_string(x = "read_count", y = "peak_count", group = "group", color = color_by)) +
      geom_path() +
      geom_point(data = lab_dt) +
      geom_text(data = lab_dt, aes(label = sample, x = read_count_lab), hjust = 0) +
      scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
      scale_y_continuous(labels = function(x)x/1e3) +
      labs(x = "reads (M)", y = "peaks (k)") +
      cowplot::theme_cowplot() +
      expand_limits(y = 0)
  }
  p_lines
}

#' plot_peak_saturation_lines.facetted
#'
#' @param cnt_dt data.table returned by load_counts()
#' @param color_by variable in cnt_dt to map to color.
#'
#' @return ggplot of peak count vs read count. facetted by sample.
#' @export
#' @import ggplot2
#' @examples
plot_peak_saturation_lines.facetted = function(cnt_dt, color_by = NULL){
  cnt_dt$group = paste(cnt_dt$out_dir, cnt_dt$sample, cnt_dt$peak_stat)
  lab_dt = data.table::rbindlist(lapply(split(cnt_dt, cnt_dt$group), function(x){
    x[x$read_count == max(x$read_count),]
  }))
  lab_dt$read_count_lab = lab_dt$read_count + 8e5
  if(is.null(color_by)){
    p = ggplot(cnt_dt, aes(x = read_count, y = peak_count, group = group)) +
      cowplot::theme_cowplot() + theme(strip.text = element_text(size = 8)) +
      scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
      scale_y_continuous(labels = function(x)x/1e3) +
      labs(x = "reads (M)", y = "peaks (k)")
    cnt_dtl = split(cnt_dt, cnt_dt$group)
    for(dt in cnt_dtl){
      p = p + annotate("path", x = dt$read_count, y = dt$peak_count, color = "gray80")
    }

    p = p +
      geom_path() +
      geom_point(data = lab_dt) +
      facet_wrap(~sample) +
      expand_limits(y = 0)
    p_lines2 = p
  }else{
    if(!color_by %in% colnames(cnt_dt)) stop(color_by, " not in colnames of cnt_dt.")
    p = ggplot(cnt_dt, aes_string(x = "read_count", y = "peak_count", group = "group", color = color_by)) +
      cowplot::theme_cowplot() + theme(strip.text = element_text(size = 8)) +
      scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
      scale_y_continuous(labels = function(x)x/1e3) +
      labs(x = "reads (M)", y = "peaks (k)")
    cnt_dtl = split(cnt_dt, cnt_dt$group)
    for(dt in cnt_dtl){
      p = p + annotate("path", x = dt$read_count, y = dt$peak_count, color = "gray80")
    }

    p = p +
      geom_path() +
      geom_point(data = lab_dt) +
      facet_wrap(~sample) +
      expand_limits(y = 0)
    p_lines2 = p
  }
  p_lines2
}

