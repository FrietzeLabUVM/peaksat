stat2flag = c("qVal" = "-q", "pVal" = "-p")

#' Title
#'
#' @slot stat character.
#' @slot stat_value numeric.
#' @slot is_PE logical.
#' @slot out_dir character.
#'
#' @return
#' @export
#'
#' @examples
setClass("peaksat_config",
         representation = list(
           stat = "character",
           stat_value = "numeric",
           is_PE = "logical",
           out_dir = "character"
         ))

#' Title
#'
#' @param stat
#' @param stat_value
#' @param is_PE
#' @param out_dir
#'
#' @return
#' @export
#'
#' @examples
#' pc = peaksat_config()
#' call_submit_script(pc, "treat.bam", "ctrl.bam")
peaksat_config = function(stat = "qVal", stat_value = .01, is_PE = FALSE, out_dir = file.path(getwd(), "peak_saturation")){

  stopifnot(stat %in% c("pVal", "qVal"))
  stopifnot(stat_value > .001)
  if(!dir.exists(out_dir))
    dir.create(out_dir, showWarnings = FALSE)
  out_dir = normalizePath(out_dir)

  new("peaksat_config", stat = stat, stat_value = stat_value, is_PE = is_PE, out_dir = out_dir)
}

make_stat_arg = function(pc){
  str = paste0(pc@stat, "_", formatC(pc@stat_value*1e3, width = 3, format = "d", flag = "0"))
  out = paste0(stat2flag[pc@stat], " ", pc@stat_value, " -o ", paste0("results_", str), ifelse(pc@is_PE, " -pe", ""))
  names(out) = str
  out
}

get_submit_script = function(){
  system.file(package = "peaksat", "extdata/submit_subsample_peak.sh")
}

get_run_script = function(){
  system.file(package = "peaksat", "extdata/run_subsample_peak.sh")
}

get_macs2_path = function(){
  system("which macs2", intern = TRUE)
}

make_submit_cmd = function(){

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
call_submit_script = function(pc, treatment, control, macs2_path = get_macs2_path(), submit_script = get_submit_script()){
  stat_arg = make_stat_arg(pc)
  paste0("bash ",
         submit_script,
         " -m ",
         macs2_path,
         " -t ",
         treatment,
         " -c ",
         control,
         " -n ",
         paste(sep = "_", sub(".bam", "", basename(treatment)), names(stat_arg)),
         " ",
         stat_arg)
}

# for(i in seq_along(t0)){
#   cmd = paste0(
#     "bash submit_subsample_peak.sh -t ",
#     t0[i],
#     " -c ",
#     c0[i],
#     " -n ",
#     paste(sep = "_", sub(".bam", "", basename(t0[i])), names(stat_arg), "vsInput"),
#     " ",
#     stat_arg,
#     ifelse(grepl("IgG", c0[i]), " -pe ", ""))
#   print(cmd)
#   system(cmd)
# }
