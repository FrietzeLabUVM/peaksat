
track_job_ids = function(jids){
  PS_OPTIONS$PS_JOB_IDS = c(getOption("PS_JOB_IDS", character()), jids)
}

#' Title
#'
#' @return
#'
#' @examples
#' treat_bam = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF10A_H3K4AC_R1/MCF10A_H3K4AC_R1.bam"
#' ctrl_bam = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF10A_input_R1/MCF10A_input_R1.bam"
#' psc = peaksat_config()
#' cmd = make_submit_cmd(psc, treat_bam, ctrl_bam)
#' sub_out = system(cmd, intern = TRUE)
#' sub_out
make_submit_cmd = function(psc,
                           treatment,
                           control,
                           submit_script = get_submit_script(),
                           hold_jid = NULL){
  stat_arg = make_stat_arg(psc, out_dir = psc@out_dir)

  if(is.null(hold_jid)){
    hold_arg = character()
  }else{
    hold_arg = paste0(" -h ", hold_jid)
  }

  paste0("bash ",
         submit_script,
         " -m ",
         psc@macs2_path,
         " -st ",
         psc@samtools_path,
         " -t ",
         treatment,
         ifelse(is.null(control),
                "",
                paste0(" -c ",
                       control)),
         ifelse(psc@is_PE,
                " -pe ",
                ""
                ),
         " -n ",
         paste(sep = "_", sub(".bam", "", basename(treatment)), names(stat_arg)),
         " ",
         stat_arg,
         " -js ",
         psc@job_scheduler,
         hold_arg)
}


#' Title
#'
#' @param psc A valid peaksat_config object
#' @param treat_bams Vector of ChIPseq bam files
#' @param ctrl_bams Vector of input bam files parallel to ChIPseq. If NULL, macs2 will be called without using an input file.
#' @param macs2_path
#' @param submit_script
#' @param await_completion
#'
#' @return
#' @export
#'
#' @examples
submit_peaksat_jobs = function(psc,
                               treat_bams,
                               ctrl_bams,
                               hold_jid_map = NULL,
                               await_completion = TRUE){
  if(!is.null(ctrl_bams)){
    if(length(ctrl_bams) != length(treat_bams)){
      if(length(ctrl_bams) != 1){
        stop("ctrl_bams must be same length as treat_bams or 1.")
      }else{
        ctrl_bams = rep(ctrl_bams, length(treat_bams))
      }
    }
  }
  if(any(!file.exists(treat_bams))){
    if(is.null(hold_jid_map)){
      stop("not all treat_bams exist!: ", paste(treat_bams[!file.exists(treat_bams)], collapse = ", "))
    }else{
      message("not all treat_bams exists, hopefully they will by time held jobs (hold_jid_map) finishes.")
      message(paste(treat_bams[!file.exists(treat_bams)], collapse = ", "))
    }

  }

  treat_bams = suppressWarnings(normalizePath(treat_bams))

  if(!is.null(ctrl_bams)){
    if(any(!file.exists(ctrl_bams))){
      if(is.null(hold_jid_map)){
        stop("not all ctrl_bams exist!: ", paste(ctrl_bams[!file.exists(ctrl_bams)], collapse = ", "))
      }else{
        message("not all ctrl_bams exists, hopefully they will by time held jobs (hold_jid_map) finishes.")
        message(paste(ctrl_bams[!file.exists(ctrl_bams)], collapse = ", "))
      }
    }
    ctrl_bams = suppressWarnings(normalizePath(ctrl_bams))
  }


  if(!is.null(hold_jid_map)){
    names(hold_jid_map) = suppressWarnings(normalizePath(names(hold_jid_map)))
  }

  cmds = sapply(seq_along(treat_bams), function(i){
    t_bam = treat_bams[i]
    c_bam = ctrl_bams[i]
    if(is.null(hold_jid_map)){
      hold_jid = NULL
    }else{
      hold_jid = hold_jid_map[t_bam]
    }
    make_submit_cmd(psc = psc, treatment = t_bam, control = c_bam, hold_jid = hold_jid)
  })
  cmd_outs = sapply(cmds, function(cmd_sub){
    qsub_str = system(cmd_sub, intern = TRUE)
    qsub_str
  })
  jids = sapply(cmd_outs, function(cmd_out){
    jid = capture_jid(cmd_out, job_scheduler = psc@job_scheduler)
    jid
  })
  PS_OPTIONS$PS_JOB_IDS = c(getOption("PS_JOB_IDS", character()), jids)

  if(await_completion){
    message(length(jids), " subset-callpeak jobs have been submitted to scheduler.\nWill await their completion before continuing.\nThis behavior can be controlled with 'await_completion = TRUE/FALSE'.\n'load_counts()' will not work properly until all jobs have finished.")
    watch_jids(jids, psc@job_scheduler)
  }else{
    message(length(jids), " subset-callpeak jobs have been submitted to scheduler.\nUse watch_jids(PS_OPTIONS$PS_JOB_IDS) to monitor.\n'load_counts()' will not work properly until all jobs have finished.")
  }

  invisible(jids)
}
