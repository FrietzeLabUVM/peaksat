#' load_counts
#'
#' @param psc peaksat_config object
#'
#' @return data.table of peak count vs read count.
#' @export
#'
#' @examples
load_counts = function(psc, min_signalValue = 1){
  wds = dir(get_result_dir(psc), full.names = TRUE)
  wds = wds[basename(wds) != "sub_logs"]
  load_counts.wd(wds, min_signalValue = min_signalValue)
}

.load_peak_counts = function(peak_count_file, min_signalValue = min_signalValue){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")

  count_files = sapply(min_signalValue, function(min_sig){
    odir = file.path(dirname(peak_count_file), paste0("counts_signalValue_", min_sig))
    dir.create(odir, showWarnings = FALSE)
    ofile = file.path(odir, basename(peak_count_file))
    ofile
  })
  if(any(!file.exists(count_files))){
    peak_file = sub(".bam.peak_count", "_peaks.narrowPeak", peak_count_file)
    if(!file.exists(peak_file)) stop("couldn't find peak file: ", peak_file)
    .peak_gr = seqsetvis::easyLoad_narrowPeak(peak_file)[[1]]
  }
  count_res = list()
  for(i in seq_along(min_signalValue)){
    min_sig = min_signalValue[i]
    odir = file.path(dirname(peak_count_file), paste0("counts_signalValue_", min_sig))
    dir.create(odir, showWarnings = FALSE)
    ofile = file.path(odir, basename(peak_count_file))
    if(file.exists(ofile)){
      dt = data.table::fread(ofile, col.names = c("peak_count", cn))
      dt$signal_cutoff = min_sig
      dt = dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
    }else{
      n_peaks = length(subset(.peak_gr, signalValue >= min_sig))
      tmp_dt = fread(peak_count_file, col.names = c("peak_count", cn))
      out_dt = copy(tmp_dt)
      out_dt$peak_count = n_peaks
      fwrite(out_dt, ofile, sep = " ")
      out_dt$signal_cutoff = min_sig
      dt = out_dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
    }
    count_res[[i]] = dt
  }
  out_dt = rbindlist(count_res)

  # for(min_sig in min_signalValue){
  #   odir = file.path(dirname(peak_count_file), paste0("counts_signalValue_", min_sig))
  #   dir.create(odir, showWarnings = FALSE)
  #   ofile = file.path(odir, basename(peak_count_file))
  #   if(file.exists(ofile)){
  #       dt = data.table::fread(ofile, col.names = c("peak_count", cn))
  #       dt$signal_cutoff = min_sig
  #       dt = dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
  #   }else{
  #     peak_file = sub(".bam.peak_count", "_peaks.narrowPeak", peak_count_file)
  #     if(!file.exists(peak_file)) stop("couldn't find peak file: ", peak_file)
  #     peak_gr = seqsetvis::easyLoad_narrowPeak(peak_file)[[1]]
  #     n_peaks = length(subset(peak_gr, signalValue >= min_sig))
  #     tmp_dt = fread(peak_count_file, col.names = c("peak_count", cn))
  #     out_dt = copy(tmp_dt)
  #     out_dt$peak_count = n_peaks
  #     fwrite(out_dt, ofile, sep = " ")
  #     out_dt$signal_cutoff = min_sig
  #     dt = out_dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
  #   }
  # }
  # if(file.exists(sub(".peak_count", ".peak_cutoff_counts", peak_count_file))){
  #   peak_cutoff_file = sub(".peak_count", ".peak_cutoff_counts", peak_count_file)
  #   dt = data.table::fread(peak_cutoff_file, col.names = c("signal_cutoff", "peak_count", cn))
  # }else{
  #   dt = data.table::fread(peak_count_file, col.names = c("peak_count", cn))
  #   dt$signal_cutoff = 1
  #   dt = dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
  # }
  out_dt[]
}

.load_counts = function(wd, min_signalValue = 1){
  read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
  peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)

  names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
  names(peak_count_files) = sub(".bam.+", "", basename(peak_count_files))

  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  rc_dt =  data.table::rbindlist(lapply(read_count_files, data.table::fread, col.names = c("read_count", cn)), idcol = "wd")
  pc_dt =  data.table::rbindlist(lapply(peak_count_files, .load_peak_counts, min_signalValue = min_signalValue), idcol = "wd")

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
load_counts.wd = function(wds, min_signalValue = 1){
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
  cnt_dt = data.table::rbindlist(lapply(wds[keep], .load_counts, min_signalValue = min_signalValue))

  cnt_dt
}
