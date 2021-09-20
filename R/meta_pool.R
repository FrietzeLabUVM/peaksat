#' make_pool_cmd
#'
#' @param pc
#' @param bams
#' @param name
#' @param pool_script
#'
#' @return
#' @export
#'
#' @examples
make_meta_pool_cmd = function(pc,
                              bams,
                              name,
                              pool_script = get_pool_script()){

  #TODO no hardcode qsub
  paste0(get_submit_command(pc@job_scheduler),
         " ",
         pool_script,
         " -st ",
         pc@samtools_path,
         " -b ",
         paste(bams, collapse = ','),
         " -n ",
         name,
         " -wd ",
         file.path(pc@out_dir, "pooled_bams")
  )
}

#' Title
#'
#' @param pc
#' @param bams
#' @param name
#' @param pool_script
#' @param skip_peaksat
#'
#' @return
#' @export
#'
#' @examples
submit_meta_pool_jobs = function(pc, bam_groups, bam_group_names = names(bams), pool_script = get_pool_script(), skip_peaksat = FALSE){
  if(!is.list(bam_groups)){
    bam_groups = list(bam_groups)
  }
  if(is.null(names(bam_groups))){
    names(bam_groups) = bam_group_names
  }
  if(is.null(names(bam_groups))){
    stop("bam_groups must be named or names must be supplied as bam_group_names")
  }
  if(any(duplicated(names(bam_groups)))){
    stop("bam_group_names must be unique")
  }
  cmds = sapply(names(bams), function(nam){
    bams = bams[[nam]]
    make_meta_pool_cmd(pc, bams, paste0(nam, "_meta"))
  })
  cmd_outs = sapply(cmds, system, intern = TRUE)
  jids = sapply(cmd_outs, function(cmd_out){
    jid = capture_jid(cmd_out)
    jid
  })
  PS_OPTIONS$PS_JOB_IDS = c(getOption("PS_JOB_IDS", character()), jids)
}
