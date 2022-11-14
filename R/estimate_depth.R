
#' Title
#'
#' @param cnt_dt
#' @param target_peaks
#' @param min_read_count
#' @param max_read_count
#'
#' @return
#' @export
#'
#' @examples
#' psc = peaksat_config("~/R/peaksat_paper/peaksat_outputs/peak_saturation_10A_H4ac_seq1_only")
#' cnt_dt = load_counts(psc)
#' est_res = estimate_depth.linear(cnt_dt[grepl("H4K8AC", sample)], target_peaks = 26e3, min_read_count = 10e6)
#' est_res$plots
#' est_res$estimates
#' est_res$estimates[, .(sample, target_Mreads = saturation_read_count/1e6, current_Mreads = current_read_count/1e6, needed_Mreads = saturation_read_count/1e6-current_read_count/1e6)]
estimate_depth.linear = function(cnt_dt, target_peaks = NULL, min_read_count = 5e6, max_read_count = 100e6){
  if(is.null(target_peaks)){
    target_peaks = .9*max(cnt_dt$peak_count)
  }
  cnt_dt[, cutoff := target_peaks]

  # saturation_dt = cnt_dt[peak_count >= cutoff, .(saturation = min(read_count), read_count, peak_count)]
  # saturation_dt = saturation_dt[saturation == read_count]
  # saturation_dt[, saturation := saturation / 1e6]
  # saturation_dt[, slope := peak_count / read_count]
  # saturation_dt[, read_count * slope]

  cnt_dt[, under_max_reads := read_count <= max_read_count]
  cnt_dt[, is_max := read_count == max(read_count[under_max_reads]), .(sample)]
  cnt_dt[, over_cutoff := peak_count > cutoff]

  min_dt = cnt_dt
  min_dt[, over_5m := read_count >= min_read_count]
  min_dt = min_dt[over_5m == TRUE]
  min_dt[, is_min := read_count == min(read_count), .(sample)]
  min_dt = min_dt[is_min == TRUE]
  min_dt[, .(sample, read_count, peak_count)]

  max_dt = cnt_dt[is_max == TRUE | over_cutoff == TRUE]
  max_dt[, selected := read_count == min(read_count), .(sample) ]
  max_dt = max_dt[selected == TRUE]
  max_dt$saturation_peak_count = target_peaks
  max_dt = merge(max_dt, min_dt[, .(zero_read_count = read_count, zero_peak_count = peak_count, sample)], by = c("sample"), all.x = TRUE)
  max_dt[is.na(zero_read_count), zero_read_count := 0]
  max_dt[is.na(zero_peak_count), zero_peak_count := 0]
  max_dt[, slope := (peak_count - zero_peak_count) / (read_count - zero_read_count)]
  max_dt[, yint := peak_count - (slope * read_count) ]
  max_dt[, saturation_read_count := (saturation_peak_count - yint) / slope]
  # max_dt[, max_read_count := max_read_count]
  # max_dt[, max_peak_count := (max_read_count - yint) / slope]
  est_dt = max_dt[, .(sample, saturation_peak_count, saturation_read_count, current_read_count = read_count, current_peak_count = peak_count, yint, slope)][order(saturation_read_count)]

  p_estimate = ggplot(cnt_dt, aes(x = read_count, y = peak_count)) +
    geom_point(data = cnt_dt, alpha = 1) +
    geom_segment(data = est_dt, lty = 1,
                 aes(x = 0,
                     y = yint,
                     xend = saturation_read_count,
                     yend = saturation_peak_count)) +
    geom_hline(yintercept = target_peaks, color = 'red') +
    geom_point(data = est_dt, aes(y = saturation_peak_count, x = saturation_read_count), color = "red") +
    geom_text(data = est_dt, aes(y = saturation_peak_count*1.07, x = saturation_read_count, label = paste0(round(saturation_read_count/1e6, 2), " (M)")), color = "red") +
    facet_wrap(~sample, scales = "free_x") +
    labs(x = "read count (M)", y = "peak count (k)") +
    scale_x_continuous(labels = function(x)x/1e6) +
    scale_y_continuous(labels = function(y)y/1e3)


  return(list(estimates = est_dt, plots = p_estimate))
}

.xy_estimate_depth.log = function(xy_dt, target_peaks, min_read_count = 5e6, max_read_count = 100e6){
  logEstimate <- lm(y~log(x-min_read_count),data=xy_dt[x >= min_read_count & x <= max_read_count])
  xvec <- seq(min_read_count+1, max_read_count,length=500)
  logpred <- predict(logEstimate, newdata=data.frame(x=xvec))

  xvec = xvec[logpred >= 0]
  logpred = logpred[logpred >= 0]

  p = ggplot(xy_dt, aes(x = x, y = y)) +
    geom_point() +
    labs(x = "read count (M)", y = "peak count (k)", title = unique(xy_dt$sample)) +
    scale_x_continuous(labels = function(x)x/1e6) +
    scale_y_continuous(labels = function(y)y/1e3) +
    geom_path(data = data.frame(x = xvec, y = logpred)) +
    geom_hline(yintercept = target_peaks, color = "red") +
    labs(subtitle = paste(round(target_peaks/1e3, 2), "k target peaks")) +
    theme(title = element_text(size = 8)) +
    expand_limits(y = 0)
  if(any(logpred >= target_peaks)){

    est_reads = xvec[min(which(logpred >= target_peaks))]
    est_peaks = logpred[min(which(logpred >= target_peaks))]
    p = p +
      annotate("point", color = "red", size = 2, x = est_reads, y = est_peaks) +
      annotate("text", color = "red", label = paste(round(est_reads/1e6, 2), "(M)"), x = est_reads, y = est_peaks*1.07, vjust = 0, hjust = .5)
  }else{
    warning("target_peaks is too high! Highest peaks estimated was ", max(logpred))
    est_reads = NA
  }
  names(est_reads) = unique(xy_dt$sample)
  return(list(estimate = est_reads, plot = p))
}

#' Estimate read depth required to reach target peaks based on logistic regression.
#'
#' @param cnt_dt data.table from load_counts
#' @param target_peaks The estimated number of peaks at saturation.
#' @param min_read_count Read counts below this value are not used when fitting regression
#' @param max_read_count Read counts above this value are not used when fitting regression
#'
#' @return A named list. Item 1 is a data.table containing estimates for saturation. Item 2 is a list of plots showing regression curves.
#' @export
#'
#' @examples
estimate_depth.log = function(cnt_dt, target_peaks = NULL, min_read_count = 5e6, max_read_count = 100e6){
  if(is.null(target_peaks)){
    target_peaks = .9*max(cnt_dt$peak_count)
  }

  xy_dtl = lapply(split(cnt_dt, cnt_dt$sample), function(dt){
    dt[, .(x = read_count, y = peak_count, sample)]
  })

  res_estimate = lapply(xy_dtl, function(xy_dt){
    # inv_log = function(x){
    #   2^x
    # }
    # expEstimate <- lm(x~(inv_log(y/2056)+min_read_count), data = xy_dt)
    # inv_log(xy_dt$y/2056)
    # predict(expEstimate, newdata = data.frame(y = seq(1, 20)*1e3))
    .xy_estimate_depth.log(xy_dt, target_peaks = target_peaks, min_read_count = min_read_count, max_read_count = max_read_count)

  })


  dt_estimate = data.table::rbindlist(lapply(res_estimate, function(x){
    estimate = x[[1]]
    xy_dt = data.table(sample = names(estimate))
    xy_dt$saturation_peak_count = target_peaks
    xy_dt$saturation_read_count = estimate
    unique(xy_dt[, .(sample, saturation_peak_count, saturation_read_count)])
  }))
  plots_estimate = lapply(res_estimate, function(x)x[[2]])
  return(list(estimates = dt_estimate, plots = plots_estimate))
}
