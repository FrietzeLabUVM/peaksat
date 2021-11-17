#' recount_peak
#'
#' @param psc
#' @param min_signalValues
#'
#' @return
#' @export
#' @import pbmcapply
#' @importFrom seqsetvis easyLoad_narrowPeak
#'
#' @examples
recount_peak = function(psc, min_signalValues = c(1, 2, 5, 10, 20, 50)){
  message("counting peaks in ", psc@out_dir)
  all_peak_files = dir(dir(dir(psc@out_dir, full.names = TRUE), full.names = TRUE), pattern = "peaks.narrowPeak$", full.names = TRUE)
  pbmcapply::pbmclapply(all_peak_files, function(peak_file){
    pc_file.orig = sub("_peaks.narrowPeak", ".bam.peak_count", peak_file)
    pc_file = sub("_peaks.narrowPeak", ".bam.peak_cutoff_counts", peak_file)
    gr = seqsetvis::easyLoad_narrowPeak(peak_file)[[1]]
    pc_dt = data.table::fread(pc_file.orig)
    pc_dt$V0 = 1
    pc_dt = pc_dt[, order(colnames(pc_dt)), with = FALSE]
    if(file.exists(pc_file))file.remove(pc_file)
    for(min_signalValue in min_signalValues){
      pc_dt$V1 = length(subset(gr, signalValue > min_signalValue))
      pc_dt$V0 = min_signalValue
      data.table::fwrite(pc_dt, file = pc_file, sep = ' ', col.names = FALSE, append = TRUE)
    }
  })
  invisible()
}
