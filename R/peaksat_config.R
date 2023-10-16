
#' Title
#'
#' @slot stat character.
#' @slot stat_value numeric.
#' @slot is_PE logical.
#' @slot out_dir character.
#' @slot macs2_path character.
#' @slot samtools_path character.
#' @slot submit_script character.
#' @slot job_scheduler character.
#' @slot noModel logical
#'
#' @slot
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
           out_dir = "character",
           macs2_path = "character",
           samtools_path = "character",
           submit_script = "character",
           job_scheduler = "character",
           noModel = "logical"
         ))

#' peaksat_config
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
#' psc = peaksat_config()
#' call_submit_script(psc, "treat.bam", "ctrl.bam")
peaksat_config = function(
  out_dir = getOption("PS_OUTDIR", file.path(getwd(), "peak_saturation")),
  stat = valid_stats$qValue,
  stat_value = .01,
  is_PE = FALSE,

  macs2_path = get_macs2_path(),
  samtools_path = get_samtools_path(),
  submit_script = get_submit_script(),
  job_scheduler = "SGE",
  noModel = FALSE
){
  if(!stat %in% valid_stats){
    stop("stat was: ", stat, ". Must be one of: ", paste(valid_stats, collapse = ", "))
  }
  stopifnot(job_scheduler %in% valid_job_schedulers)
  if(!dir.exists(out_dir))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_dir = normalizePath(out_dir)

  run_script = file.path(dirname(submit_script), "run_subsample_peak.sh")
  if(!file.exists(submit_script)){
    stop("submit_script not found!: ", submit_script)
  }
  if(!file.exists(run_script)){
    stop("run_script not found!: ", run_script)
  }

  new("peaksat_config",
      stat = stat,
      stat_value = stat_value,
      is_PE = is_PE,
      out_dir = out_dir,
      macs2_path = macs2_path,
      samtools_path = samtools_path,
      submit_script = submit_script,
      job_scheduler = job_scheduler,
      noModel = noModel)
}

#' peaksat_config.example
#'
#' @return A peaksat_config object using data packaged with peaksat.
#' @export
#'
#' @examples
#' psc = peaksat_config.example()
peaksat_config.example = function(){
  ex_dir = system.file(package = "peaksat", "extdata/example_peaksat_output")
  psc = peaksat_config(out_dir = ex_dir, stat_value = .01)
  psc
}
