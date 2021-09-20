
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
                           submit_script = get_submit_script(),
                           hold_jid = NULL){
  stat_arg = make_stat_arg(pc, out_dir = pc@out_dir)

  hold_arg = ifelse(is.null(hold_jid),
                    character(),
                    paste0(" -h ", hold_jid)
  )

  paste0("bash ",
         submit_script,
         " -m ",
         pc@macs2_path,
         " -st ",
         pc@samtools_path,
         " -t ",
         treatment,
         " -c ",
         control,
         " -n ",
         paste(sep = "_", sub(".bam", "", basename(treatment)), names(stat_arg)),
         " ",
         stat_arg,
         " -js ",
         pc@job_scheduler,
         hold_arg)
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
                               hold_jid_map = NULL,
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

  treat_bams = normalizePath(treat_bams)
  ctrl_bams = normalizePath(ctrl_bams)
  if(!is.null(hold_jid_map)){
    stopifnot(file.exists(names(hold_jid_map)))
    names(hold_jid_map) = normalizePath(names(hold_jid_map))
  }

  cmds = sapply(seq_along(treat_bams), function(i){
    t_bam = treat_bams[i]
    c_bam = ctrl_bams[i]
    hold_jid = ifelse(is.null(hold_jid_map), NULL, hold_jid_map[t_bam])
    make_submit_cmd(pc = pc, treatment = t_bam, control = c_bam, hold_jid = hold_jid)
  })
  cmd_outs = sapply(cmds, function(cmd_sub){
    qsub_str = system(cmd_sub, intern = TRUE)
    qsub_str
  })
  jids = sapply(cmd_outs, function(cmd_out){
    jid = capture_jid(cmd_out)
    jid
  })
  PS_OPTIONS$PS_JOB_IDS = c(getOption("PS_JOB_IDS", character()), jids)

  if(await_completion){
    message(length(jids), " subset-callpeak jobs have been submitted to scheduler.\nWill await their completion before continuing.\nThis behavior can be controlled with 'await_completion = TRUE/FALSE'.\n'load_counts()' will not work properly until all jobs have finished.")
    watch_jids(jids)
  }else{
    message(length(jids), " subset-callpeak jobs have been submitted to scheduler.\nUse watch_jids() to monitor.\n'load_counts()' will not work properly until all jobs have finished.")
  }

  invisible(jids)
}
