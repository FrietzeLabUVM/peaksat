load_counts(ps){
 wds = dir(ps@out_dir, full.names = TRUE)
}

#' Title
#'
#' @param wds
#'
#' @return
#' @export
#' @import data.table
#' @examples
load_counts.wd = function(wds){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")

  .load_counts = function(wd){
    read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
    peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)


    names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
    names(peak_count_files) = sub(".bam.+", "", basename(peak_count_files))

    rc_dt = lapply(read_count_files, data.table::fread, col.names = c("read_count", cn)) %>%
      data.table::rbindlist(idcol = "wd")
    pc_dt = lapply(peak_count_files, data.table::fread, col.names = c("peak_count", cn)) %>%
      data.table::rbindlist(idcol = "wd")


    cnt_dt = merge(rc_dt, pc_dt[, .(name, peak_count)], by = "name")

    cnt_dt[, c("fraction", "rep") := data.table::tstrsplit(name, "\\.", keep = 2:3)]
    cnt_dt[, fraction := as.numeric(fraction)/100]

    cnt_dt$sample = sub("peak_saturation.", "", basename(wd))

    cnt_dt[order(read_count)]
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



ggplot(cnt_dt, aes(x = read_count, y = peak_count, group = sample, color = cell)) +
  geom_path() +
  geom_point(data = lab_dt) +
  geom_text(data = lab_dt, aes(label = paste(cell, mark, rep), x = read_count + 800e3), hjust = 0, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
  scale_y_continuous(labels = function(x)x/1e3) +
  labs(x = "reads (M)", y = "peaks (k)") +
  cowplot::theme_cowplot()

plot_peak_saturation_lines = function(cnt_dt){
  p_lines = ggplot(cnt_dt, aes(x = read_count, y = peak_count, group = sample)) +
    geom_path() +
    geom_point(data = lab_dt) +
    geom_text(data = lab_dt, aes(label = sample, x = read_count + 800e3), hjust = 0) +
    # scale_x_continuous(expand = c(.1, 2e9))
    scale_x_continuous(expand = expansion(c(.1, .5)), labels = function(x)x/1e6) +
    scale_y_continuous(labels = function(x)x/1e3) +
    labs(x = "reads (M)", y = "peaks (k)") +
    cowplot::theme_cowplot()
  p_lines
}

plot_peak_saturation_lines.facetted = function(cnt_dt){
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
    facet_wrap(~sample)
  p_lines2 = p
  p_lines2
}

plot_peak_saturation_lines(cnt_dt)
plot_peak_saturation_lines.facetted(cnt_dt)
