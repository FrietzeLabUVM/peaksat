
get_str = function(psc){
  pow10 = max(3, -log10(psc@stat_value)+1)
  paste0(psc@stat, "_", formatC(psc@stat_value*10^pow10, width = pow10, format = "d", flag = "0"))
}

make_stat_arg = function(psc, out_dir = getOption("PS_OUTDIR", getwd())){
  str = get_str(psc)
  out = paste0(stat2flag[psc@stat], " ", psc@stat_value, " -o ", get_result_dir(psc), ifelse(psc@is_PE, " -pe", ""))
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
#' @param psc
#'
#' @return
#' @export
#'
#' @examples
get_result_dir = function(psc){
  str = get_str(psc)
  file.path(psc@out_dir, paste0("results_", str))
}

get_pool_dir = function(psc){
  str = get_str(psc)
  file.path(psc@out_dir, paste0("pooled_bams"))
}
