#' load_counts
#'
#' @param psc
#'
#' @return
#' @export
#'
#' @examples
load_counts = function(psc){
 wds = dir(get_result_dir(psc), full.names = TRUE)
 wds = wds[basename(wds) != "sub_logs"]
 load_counts.wd(wds)
}

#' Title
#'
#' @param wds
#'
#' @return
#' @export
#'
#' @examples
load_counts.wd = function(wds){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")

  .load_counts = function(wd){
    read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
    peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)

    names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
    names(peak_count_files) = sub(".bam.+", "", basename(peak_count_files))

    rc_dt =  data.table::rbindlist(lapply(read_count_files, data.table::fread, col.names = c("read_count", cn)), idcol = "wd")
    pc_dt =  data.table::rbindlist(lapply(peak_count_files, data.table::fread, col.names = c("peak_count", cn)), idcol = "wd")

    cnt_dt = merge(rc_dt, pc_dt[, c("name", "peak_count")], by = "name")

    data.table::set(cnt_dt, j = c("fraction", "rep"), value = data.table::tstrsplit(cnt_dt$name, "\\.", keep = 2:3))
    cnt_dt$fraction = as.numeric(cnt_dt$fraction)/100

    cnt_dt$sample = sub("peak_saturation.", "", basename(wd))

    res_dir = basename(dirname(wds))
    stat_suff = sub("results_", "", res_dir)

    samp_names = sub("peak_saturation.", "", basename(wds))

    clean_samp_names = sapply(seq_along(samp_names), function(i){
      sub(paste0("_", stat_suff[i]), "", samp_names[i])
    })

    names(stat_suff) = samp_names
    names(clean_samp_names) = samp_names

    cnt_dt[, peak_stat := stat_suff[sample]]
    cnt_dt[, sample := clean_samp_names[sample]]
    cnt_dt[, c("stat_name", "stat_value") := tstrsplit(peak_stat, "_") ]
    cnt_dt[, stat_value := as.numeric(stat_value)/10^nchar(stat_value)]
    cnt_dt[order(cnt_dt$read_count),]
  }

  wds = wds[basename(wds) != "sub_logs"]
  keep = sapply(wds, function(wd){
    read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
    peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)
    if(length(read_count_files) == 0 || length(peak_count_files) == 0){
      return(FALSE)
    }
    if(length(read_count_files) != length(peak_count_files)){
      return(FALSE)
    }
    return(TRUE)
  })
  if(any(!keep)){
    warning("Not all results appear to be complete.")
  }
  if(!any(keep)) stop ("No valid results found.")
  cnt_dt = data.table::rbindlist(lapply(wds[keep], .load_counts))

  cnt_dt
}

#' Title
#'
#' @param cnt_dt
#'
#' @return
#' @export
#' @import ggplot2
#' @examples
plot_peak_saturation_lines = function(cnt_dt){
  lab_dt = data.table::rbindlist(lapply(split(cnt_dt, cnt_dt$sample), function(x){
    x[x$read_count == max(x$read_count),]
  }))
  lab_dt$read_count_lab = lab_dt$read_count + 8e5
  p_lines = ggplot(cnt_dt, aes(x = read_count, y = peak_count, group = sample)) +
    geom_path() +
    geom_point(data = lab_dt) +
    geom_text(data = lab_dt, aes(label = sample, x = read_count_lab), hjust = 0) +
    scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
    scale_y_continuous(labels = function(x)x/1e3) +
    labs(x = "reads (M)", y = "peaks (k)") +
    cowplot::theme_cowplot() +
    expand_limits(y = 0)
  p_lines
}

#' Title
#'
#' @param cnt_dt
#'
#' @return
#' @export
#' @import ggplot2
#' @examples
plot_peak_saturation_lines.facetted = function(cnt_dt){
  lab_dt = data.table::rbindlist(lapply(split(cnt_dt, cnt_dt$sample), function(x){
    x[x$read_count == max(x$read_count),]
  }))
  lab_dt$read_count_lab = lab_dt$read_count + 8e5
  p = ggplot(cnt_dt, aes(x = read_count, y = peak_count, group = sample)) +
    cowplot::theme_cowplot() + theme(strip.text = element_text(size = 8)) +
    scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
    scale_y_continuous(labels = function(x)x/1e3) +
    labs(x = "reads (M)", y = "peaks (k)")
  cnt_dtl = split(cnt_dt, cnt_dt$sample)
  for(dt in cnt_dtl){
    p = p + annotate("path", x = dt$read_count, y = dt$peak_count, color = "gray80")
  }

  p = p +
    geom_path() +
    geom_point(data = lab_dt) +
    facet_wrap(~sample) +
    expand_limits(y = 0)
  p_lines2 = p
  p_lines2
}

