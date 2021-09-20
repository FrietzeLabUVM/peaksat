
get_str = function(pc){
  paste0(pc@stat, "_", formatC(pc@stat_value*1e3, width = 3, format = "d", flag = "0"))
}

make_stat_arg = function(pc, out_dir = getOption("PS_OUTDIR", getwd())){
  str = get_str(pc)
  out = paste0(stat2flag[pc@stat], " ", pc@stat_value, " -o ", get_result_dir(pc), ifelse(pc@is_PE, " -pe", ""))
  names(out) = str
  out
}

get_submit_script = function(){
  system.file(package = "peaksat", "extdata/submit_subsample_peak.sh")
}

get_pool_script = function(){
  system.file(package = "peaksat", "extdata/run_pool_bams.sh")
}

get_run_script = function(){
  system.file(package = "peaksat", "extdata/run_subsample_peak.sh")
}

get_macs2_path = function(){
  getOption("PS_MACS2_PATH", system("which macs2", intern = TRUE))
}

get_samtools_path = function(){
  getOption("PS_SAMTOOLS_PATH", system("which samtools", intern = TRUE))
}

#' Title
#'
#' @param pc
#'
#' @return
#' @export
#'
#' @examples
get_result_dir = function(pc){
  str = get_str(pc)
  file.path(pc@out_dir, paste0("results_", str))
}

