

#' load_overlap_convergence
#'
#' @template param_sqc
#' @param a
#' @param b
#' @param min_signalValue
#'
#' @return
#' @export
#' @rdname convergence
#' @examples
#' psc = peaksat_config.example()
#' cnv_dt = load_convergence(psc, "rep1", "rep2")
#' plot_convergence_bars(cnv_dt, "rep1", "rep2")
#' plot_convergence_trends(cnv_dt, "rep1", "rep2")
load_overlap_convergence = function(psc, a, b, min_signalValue = 1){
  all_res = dir(get_result_dir(psc), full.names = TRUE)

  avail_res = show_available_results(psc)

  a_dir = avail_res[grepl(a, avail_res)]
  b_dir = avail_res[grepl(b, avail_res)]
  if(length(min_signalValue) != 1){
    stop("Currently only a single min_signalValue cutoff is supported.")
  }
  if(length(a_dir) < 1){
    stop(a, " did match any results, use show_available_results to see valid options.")
  }
  if(length(b_dir) < 1){
    stop(b, " did match any results, use show_available_results to see valid options.")
  }
  if(length(a_dir) > 1){
    stop(a, " matched too many results:\n", paste(a_dir, collapse = "\n"))
  }
  if(length(b_dir) > 1){
    stop(b, " matched too many results:\n", paste(b_dir, collapse = "\n"))
  }

  cnt_dt = load_counts(psc, selection = c(a_dir, b_dir), min_signalValue = min_signalValue)
  cnt_dt = cnt_dt[, .(name, read_count, sample)]
  cnt_dt = dcast(cnt_dt, name~sample, value.var = "read_count")
  setnames(cnt_dt, c(a_dir, b_dir), c("A_read_count", "B_read_count"))

  peak_grs = load_peaks(psc, selection = c(a_dir, b_dir), min_signalValue = min_signalValue)
  a_grs = peak_grs[[1]][[a_dir]]
  b_grs = peak_grs[[1]][[b_dir]]
  all_olap_grs = lapply(seq_along(a_grs), function(i){
    table(seqsetvis::ssvFactorizeMembTable(seqsetvis::ssvOverlapIntervalSets(list(A = a_grs[[i]], B = b_grs[[i]])))$group)
  })
  names(all_olap_grs) = names(a_grs)
  olap_dt = rbindlist(lapply(all_olap_grs, function(x){
    dt = data.table(x)
    setnames(dt, c("group", "count"))
    dt
  }), idcol = "name")
  olap_dt = merge(olap_dt, cnt_dt, by = c("name"))

  olap_dt = olap_dt[group != "none"]

  olap_dt[,  center_width := count[group == "A & B"] , .(name)]
  olap_dt[, xmin := -center_width / 2]
  olap_dt[, xmax := center_width / 2]
  olap_dt[group == "A", xmax := xmin]
  olap_dt[group == "A", xmin := xmin - count]
  olap_dt[group == "B", xmin := xmax]
  olap_dt[group == "B", xmax := xmax + count]

  olap_dt[, fraction_str := tstrsplit(name, "\\.", keep = 2)]
  olap_dt[, fraction := as.numeric(fraction_str)/100]

  olap_dt[, ymax := fraction + .05]
  olap_dt[, ymin := fraction - .05]

  olap_dt[, total := sum(count), .(name)]
  olap_dt[, count_fraction := count / total]
  olap_dt[]
}

#' plot_overlap_convergence_bars
#'
#' @param olap_dt output from load_overlap_convergence function
#' @rdname convergence
#'
#' @return ggplot showing stacked bars of overlaps at different subsample levels.
#' @export
plot_overlap_convergence_bars = function(olap_dt, a, b){
  ggplot(olap_dt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group)) +
    geom_rect(color = "black") +
    labs(
      y = "fraction of reads",
      x = "number of peaks",
      title = paste(sep = "\n",
                    paste("A =", a),
                    paste("B =", b)
      )
    )
}



#' plot_overlap_convergence_trends
#'
#' @param olap_dt
#' @rdname convergence
#'
#' @return ggplot showing line plots of overlap trends at different subsample levels.
#' @export
plot_overlap_convergence_trends = function(olap_dt, a, b){
  ggplot(olap_dt, aes(x = fraction, y = count_fraction, color = group))  +
    geom_path() +
    geom_point() +
    # geom_point(aes(size = count_fraction)) +
    labs(y = "fraction of peaks", x = "subset of reads",
         title = paste(sep = "\n",
                       paste("A =", a),
                       paste("B =", b)
         ))
}
