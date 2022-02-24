#' load_counts
#'
#' @param psc peaksat_config object
#' @param min_signalValue
#'
#' @return data.table of peak count vs read count.
#' @export
#'
#' @examples
load_counts = function(psc, min_signalValue = 1, selection = NULL){
  .general_load_info(psc,
                     min_signalValue = min_signalValue,
                     type = "count",
                     selection = selection)
}

#' load_widths
#'
#' @param psc peaksat_config object
#' @param min_signalValue
#'
#' @return data.table of peak count vs read count.
#' @export
#'
#' @examples
load_widths = function(psc, min_signalValue = 1, selection = NULL){
  .general_load_info(psc,
                     min_signalValue = min_signalValue,
                     type = "width",
                     selection = selection)
}

#' load_novelty
#'
#' @param psc peaksat_config object
#' @param min_signalValue
#'
#' @return data.table of peak count vs read count.
#' @export
#'
#' @examples
load_novelty = function(psc, min_signalValue = 1, selection = NULL){
  .general_load_info(psc,
                     min_signalValue = min_signalValue,
                     type = "novelty",
                     selection = selection)
}

#' load_peaks
#'
#' @param psc peaksat_config object
#' @param min_signalValue
#'
#' @return list of GRanges of peaks filtered by min_signalValue.
#' @export
#'
#' @examples
load_peaks = function(psc, min_signalValue = 1, selection = NULL){
  .general_load_info(psc,
                     min_signalValue = min_signalValue,
                     type = "peaks",
                     selection = selection)
}

.general_load_info = function(psc, min_signalValue = 1, type = c("count", "width", "novelty", "peaks")[1], selection = NULL){
  if(is.null(selection)){
    wds = dir(get_result_dir(psc), full.names = TRUE)
    wds = wds[basename(wds) != "sub_logs"]
  }else{
    res_dir = get_result_dir(psc)
    stat_str = sub("results_", "", basename(res_dir))
    wds = file.path(res_dir, paste0("peak_saturation.", selection, "_", stat_str))
    if(!all(dir.exists(wds))) stop("Invalid result directory in selection:\n", paste(wds[!all(dir.exists(wds))], collapse = "\n"))
  }

  .general_load_info.wd(wds, min_signalValue = min_signalValue, type = type)
}

#' .general_load_info.wd
#'
#' @param wds directories to load count data from.
#'
#' @return data.table of peak count vs read count.
#'
#' @examples
.general_load_info.wd = function(wds, min_signalValue = 1, type = c("count", "width", "novelty", "peaks")[1]){
  wds = wds[basename(wds) != "sub_logs"]
  keep = sapply(wds, function(wd){
    read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
    np_files = dir(wd, pattern = "narrowPeak$", full.names = TRUE)
    if(length(read_count_files) == 0 || length(np_files) == 0){
      return(FALSE)
    }
    if(length(read_count_files) != length(np_files)){
      return(FALSE)
    }
    return(TRUE)
  })
  if(any(!keep)){
    warning("Not all results appear to be complete.")
  }
  if(!any(keep)) stop ("No valid results found.")


  peak_FUN = switch(type,
                    count = {
                      .load_peak_counts
                    },
                    width = {
                      .load_peak_widths
                    },
                    novelty = {
                      .load_peak_novelty
                    },
                    peaks = {
                      .load_peak_peaks
                    }, {
                      NULL
                    })
  if(is.null(peak_FUN)){
    stop("Unrecognized type for peak_FUN: ", type)
  }

  if(type == "peaks"){
    all_peaks = pbapply::pblapply(wds[keep], function(wd){
      np_files = dir(wd, pattern = "narrowPeak$", full.names = TRUE)
      names(np_files) = sub("_peaks.narrowPeak", "", basename(np_files))
      seqsetvis::easyLoad_narrowPeak(np_files)
    })
    names(all_peaks) = sub(paste0("_", names(make_stat_arg(psc))), "", sub("^peak_saturation.", "", basename(wds[keep])))
    #need to invert nested list, currently sample[sig_value], need sig_value[sample]
    names(min_signalValue) = paste0("sigValue_", min_signalValue)
    out_list = lapply(min_signalValue, function(min_sig){
      lapply(all_peaks, function(sample_peaks){
        lapply(sample_peaks, function(peak_gr){
          subset(peak_gr, signalValue >= min_sig)
        })
      })
    })
    return(out_list)
  }

  if(.Platform$OS.type == "windows" || getOption("mc.cores", 1) == 1) {
    cnt_dtl = pbapply::pblapply(wds[keep], .assemble_peak_info, min_signalValue = min_signalValue, info_name = type, peak_FUN = peak_FUN)

  } else {
    cnt_dtl = pbmcapply::pbmclapply(wds[keep], .assemble_peak_info, min_signalValue = min_signalValue, info_name = type, peak_FUN = peak_FUN)
  }



  cnt_dt = data.table::rbindlist(cnt_dtl)

  cnt_dt[]
}

.get_template = function(np_file){
  suffix = "count"
  odir = dirname(np_file)
  ofile = file.path(odir, sub("_peaks.narrowPeak", paste0(".bam.peak_", suffix), basename(np_file)))
  ofile
}

.get_ofile = function(np_file, min_sig, suffix = c("count", "width", "novelty")[1]){
  odir = file.path(dirname(np_file), paste0("counts_signalValue_", min_sig))
  dir.create(odir, showWarnings = FALSE)
  ofile = file.path(odir, sub("_peaks.narrowPeak", paste0(".bam.peak_", suffix), basename(np_file)))
  ofile
}

.assemble_peak_info = function(wd, min_signalValue = 1, info_name = c("count", "width", "novelty")[1], peak_FUN = .load_peak_counts){
  read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
  np_files = dir(wd, pattern = "narrowPeak$", full.names = TRUE)

  names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
  names(np_files) = sub("_peaks.narrowPeak", "", basename(np_files))

  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  rc_dt =  data.table::rbindlist(lapply(read_count_files, data.table::fread, col.names = c("read_count", cn)), idcol = "wd")
  if(info_name == "novelty"){
    ofile = .get_ofile(np_files[1], min_sig = min_sig, suffix = "novelty")
    if(file.exists(ofile)){
      pc_dt = fread(ofile)
    }else{
      pc_grs = lapply(np_files, peak_FUN, min_signalValue = min_signalValue)

      temp_files = sapply(np_files, .get_template)
      tmp_dt = rbindlist(lapply(temp_files, fread, col.names = c("V1", cn)))
      tmp_dt$V1 = NULL
      out_dt = copy(tmp_dt)
      out_dt$peak_novelty = -1

      peak_groups = lapply(seq_along(pc_grs[[1]]), function(i){
        lapply(pc_grs, function(x){
          x[[i]]
        })
      })
      names(peak_groups) = names(pc_grs[[1]])
      for(nam in names(peak_groups)){
        min_sig = as.numeric(nam)
        grp = peak_groups[[nam]]
        o = order(as.numeric(sapply(strsplit(names(grp), "\\."), function(x)x[2])))
        grp = grp[o]

        so_far = GRanges()
        for(i in seq_along(grp)){
          novel_peaks = subsetByOverlaps(grp[[i]], so_far, invert = TRUE)
          n_peaks = length(novel_peaks)
          message(n_peaks)
          out_dt$peak_novelty[i] = n_peaks
          so_far = suppressWarnings({
            reduce(c(so_far, grp[[i]]))
          })
        }
        out_dt$signal_cutoff = min_sig
        pc_dt = out_dt[, c("signal_cutoff", "peak_novelty", cn), with = FALSE]
        fwrite(pc_dt, ofile, sep = " ")
      }
      pc_dt
    }

  }else{
    pc_dt =  data.table::rbindlist(lapply(np_files, peak_FUN, min_signalValue = min_signalValue), idcol = "wd")
  }


  cnt_dt = merge(rc_dt, pc_dt[, c("name", paste0("peak_", info_name), "signal_cutoff"), with = FALSE], by = "name")

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

#loads a vector of peak_count_files
.load_peak_peaks = function(np_file, min_signalValue = 1){
  message("Loading peak file: ", np_file)
  if(!file.exists(np_file)) stop("couldn't find peak file: ", np_file)
  .peak_gr = seqsetvis::easyLoad_narrowPeak(np_file)[[1]]
  peak_grs = list()
  for(i in seq_along(min_signalValue)){
    min_sig = min_signalValue[i]
    peak_grs[[i]] = subset(.peak_gr, signalValue >= min_sig)
  }
  peak_grs
}

#loads a vector of peak_count_files
.load_peak_counts = function(np_file, min_signalValue = 1){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  count_files = sapply(min_signalValue, function(min_sig){
    .get_ofile(np_file, min_sig, "count")
  })
  if(any(!file.exists(count_files))){
    message("Loading peak file for counting: ", np_file)
    if(!file.exists(np_file)) stop("couldn't find peak file: ", np_file)
    .peak_gr = seqsetvis::easyLoad_narrowPeak(np_file)[[1]]
  }else{
    message("Using saved peak counts for: ", np_file)
  }
  count_res = list()
  for(i in seq_along(min_signalValue)){
    min_sig = min_signalValue[i]
    ofile = .get_ofile(np_file, min_sig, "count")
    if(file.exists(ofile)){
      dt = data.table::fread(ofile, col.names = c("peak_count", cn))
      dt$signal_cutoff = min_sig
      dt = dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
    }else{
      n_peaks = length(subset(.peak_gr, signalValue >= min_sig))
      tmp_dt = fread(.get_template(np_file), col.names = c("peak_count", cn))
      out_dt = copy(tmp_dt)
      out_dt$peak_count = n_peaks
      fwrite(out_dt, ofile, sep = " ")
      out_dt$signal_cutoff = min_sig
      dt = out_dt[, c("signal_cutoff", "peak_count", cn), with = FALSE]
    }
    count_res[[i]] = dt
  }
  out_dt = rbindlist(count_res)

  out_dt[]
}

#loads a vector of peak_width_files
.load_peak_widths = function(np_file, min_signalValue = 1){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  width_files = sapply(min_signalValue, function(min_sig){
    .get_ofile(np_file, min_sig, "width")
  })
  if(any(!file.exists(width_files))){
    message("Loading peak file for totaling width: ", np_file)
    if(!file.exists(np_file)) stop("couldn't find peak file: ", np_file)
    .peak_gr = seqsetvis::easyLoad_narrowPeak(np_file)[[1]]
  }else{
    message("Using saved peak widths for: ", np_file)
  }
  width_res = list()
  for(i in seq_along(min_signalValue)){
    min_sig = min_signalValue[i]
    ofile = .get_ofile(np_file, min_sig, "width")
    if(file.exists(ofile)){
      dt = data.table::fread(ofile, col.names = c("peak_width", cn))
      dt$signal_cutoff = min_sig
      dt = dt[, c("signal_cutoff", "peak_width", cn), with = FALSE]
    }else{
      peak_width = sum(width(subset(.peak_gr, signalValue >= min_sig)))
      tmp_dt = fread(.get_template(np_file), col.names = c("peak_width", cn))
      out_dt = copy(tmp_dt)
      out_dt$peak_width = peak_width
      fwrite(out_dt, ofile, sep = " ")
      out_dt$signal_cutoff = min_sig
      dt = out_dt[, c("signal_cutoff", "peak_width", cn), with = FALSE]
    }
    width_res[[i]] = dt
  }
  out_dt = rbindlist(width_res)

  out_dt[]
}

.load_peak_novelty = function(np_file, min_signalValue = 1){
  message("Loading peak file for novelty calculation: ", np_file)
  if(!file.exists(np_file)) stop("couldn't find peak file: ", np_file)
  .peak_gr = seqsetvis::easyLoad_narrowPeak(np_file)[[1]]
  peak_grs = list()
  for(i in seq_along(min_signalValue)){
    min_sig = min_signalValue[i]
    peak_grs[[as.character(min_sig)]] = subset(.peak_gr, signalValue >= min_sig)
  }
  peak_grs
}

#' show_available_results
#'
#' @param psc
#'
#' @return
#' @export
#'
#' @examples
show_available_results = function(psc){
  res = dir(get_result_dir(psc))
  res = res[res != "sub_logs"]
  res = sub("^peak_saturation.", "", res)
  res = sub(sub("results", "", basename(get_result_dir(psc))), "", res)
  res
}
