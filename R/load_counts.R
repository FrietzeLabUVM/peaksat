#' load_counts
#'
#' @param psc peaksat_config object
#'
#' @return data.table of peak count vs read count.
#' @export
#'
#' @examples
load_counts = function(psc){
  wds = dir(get_result_dir(psc), full.names = TRUE)
  wds = wds[basename(wds) != "sub_logs"]
  load_counts.wd(wds)
}

peak_count_file = "/slipstream/home/joeboyd/R/peaksat_paper/peak_saturation/results_qValue_010/peak_saturation.CD34-01517_H3K27me3_R1_qValue_010/subset.010.2.bam.peak_count"
.load_peak_counts = function(peak_count_file){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  if(file.exists(sub(".peak_count", ".peak_cutoff_counts", peak_count_file))){
    peak_cutoff_file = sub(".peak_count", ".peak_cutoff_counts", peak_count_file)
    dt = data.table::fread(peak_cutoff_file, col.names = c("signal_cutoff", "peak_count", cn))
  }else{
    dt = data.table::fread(peak_count_file, col.names = c("peak_count", cn))
    dt$signal_cutoff = 1
    dt = dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
  }
  dt
}

.load_counts = function(wd){
  read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
  peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)

  names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
  names(peak_count_files) = sub(".bam.+", "", basename(peak_count_files))

  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  rc_dt =  data.table::rbindlist(lapply(read_count_files, data.table::fread, col.names = c("read_count", cn)), idcol = "wd")
  pc_dt =  data.table::rbindlist(lapply(peak_count_files, .load_peak_counts), idcol = "wd")

  cnt_dt = merge(rc_dt, pc_dt[, c("name", "peak_count", "signal_cutoff")], by = "name")

  data.table::set(cnt_dt, j = c("fraction", "rep"), value = data.table::tstrsplit(cnt_dt$name, "\\.", keep = 2:3))
  cnt_dt$fraction = as.numeric(cnt_dt$fraction)/100

  cnt_dt$sample = sub("peak_saturation.", "", basename(wd))

  res_dir = basename(dirname(wd))
  stat_suff = sub("results_", "", res_dir)

  samp_names = sub("peak_saturation.", "", basename(wd))

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

#' load_counts.wd
#'
#' @param wds directories to load count data from.
#'
#' @return data.table of peak count vs read count.
#' @export
#'
#' @examples
load_counts.wd = function(wds){
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
