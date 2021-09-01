

#' Title
#'
#' @slot stat character.
#' @slot stat_value numeric.
#' @slot is_PE logical.
#' @slot out_dir character.
#' @slot macs2_path character
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
           submit_script = "character"
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
peaksat_config = function(stat = valid_stats$qValue,
                          stat_value = .01,
                          is_PE = FALSE,
                          out_dir = getOption("PS_OUTDIR", file.path(getwd(), "peak_saturation")),
                          macs2_path = get_macs2_path(),
                          submit_script = get_submit_script()){

  stopifnot(stat %in% valid_stats)
  stopifnot(stat_value > .001)
  if(!dir.exists(out_dir))
    dir.create(out_dir, showWarnings = FALSE)
  out_dir = normalizePath(out_dir)

  run_script = file.path(dirname(submit_script), "run_subsample_peak.sh")
  if(!file.exists(submit_script)){
    stop("submit_script not found!: ", submit_script)
  }
  if(!file.exists(run_script)){
    stop("run_script not found!: ", run_script)
  }

  new("peaksat_config", stat = stat, stat_value = stat_value, is_PE = is_PE, out_dir = out_dir, macs2_path = macs2_path, submit_script = submit_script)
}

make_stat_arg = function(pc, out_dir = getOption("PS_OUTDIR", getwd())){
  str = paste0(pc@stat, "_", formatC(pc@stat_value*1e3, width = 3, format = "d", flag = "0"))
  out = paste0(stat2flag[pc@stat], " ", pc@stat_value, " -o ", file.path(paste0("results_", str)), ifelse(pc@is_PE, " -pe", ""))
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
  getOption("PS_MACS2_PATH", system("which macs2", intern = TRUE))
}

#' Title
#'
#' @return
#'
#' @examples
#' treat_bam = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF10A_H3K4AC_R1/MCF10A_H3K4AC_R1.bam"
#' ctrl_bam = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF10A_input_R1/MCF10A_input_R1.bam"
#' pc = peaksat_config()
#' cmd = make_submit_cmd(pc, treat_bam, ctrl_bam)
#' sub_out = system(cmd, intern = TRUE)
#' sub_out
make_submit_cmd = function(pc,
                           treatment,
                           control,
                           submit_script = get_submit_script()){
  stat_arg = make_stat_arg(pc, out_dir = pc@out_dir)
  paste0("bash ",
         submit_script,
         " -m ",
         pc@macs2_path,
         " -t ",
         treatment,
         " -c ",
         control,
         " -n ",
         paste(sep = "_", sub(".bam", "", basename(treatment)), names(stat_arg)),
         " ",
         stat_arg)
}

#' Title
#'
#' @param pc
#' @param treat_bams
#' @param ctrl_bams
#' @param macs2_path
#' @param submit_script
#' @param await_completion
#'
#' @return
#' @export
#'
#' @examples
submit_peaksat_jobs = function(pc,
                               treat_bams,
                               ctrl_bams,
                               await_completion = TRUE){
  if(length(ctrl_bams) != length(treat_bams)){
    if(length(ctrl_bams) != 1){
      stop("ctrl_bams must be same length as treat_bams or 1.")
    }else{
      ctrl_bams = rep(ctrl_bams, length(treat_bams))
    }
  }
  if(any(!file.exists(treat_bams))){
    stop("not all treat_bams exist!: ", paste(treat_bams[!file.exists(treat_bams)], collapse = ", "))
  }
  if(any(!file.exists(ctrl_bams))){
    stop("not all treat_bams exist!: ", paste(ctrl_bams[!file.exists(ctrl_bams)], collapse = ", "))
  }
  cmds = sapply(seq_along(treat_bams), function(i){
    t_bam = treat_bams[i]
    c_bam = ctrl_bams[i]
    make_submit_cmd(pc = pc, treatment = t_bam, control = c_bam)
  })
  if(await_completion){
    jids = qsub_and_wait(cmds)
  }else{
    jids = sapply(cmds, function(cmd_sub){
      qsub_str = system(cmd_sub, intern = TRUE)
      jid = capture_jid(qsub_str)
      jid
    })
    PS_OPTIONS$PS_JOB_IDS = c(getOption("PS_JOB_IDS", character()), jids)
  }
  invisible(jids)
}
