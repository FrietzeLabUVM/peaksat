#' make_pool_cmd
#'
#' @param psc
#' @param bams
#' @param name
#' @param pool_script
#'
#' @return
#' @export
#'
#' @examples
make_meta_pool_cmd = function(psc,
                              bams,
                              name,
                              pool_script = get_pool_script()){
  p_bam_dir = file.path(get_pool_dir(psc), "pooled_bams")
  log_dir = file.path(p_bam_dir, "sub_logs")
  dir.create(p_bam_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  paste0(get_submit_command(psc@job_scheduler),
         " -o ", log_dir,
         " -e ", log_dir,
         " ",
         pool_script,
         " -st ",
         psc@samtools_path,
         " -b ",
         paste(bams, collapse = ','),
         " -n ",
         name,
         " -wd ",
         file.path(psc@out_dir, "pooled_bams")
  )
}

get_meta_bam_file = function(psc, bam_group_name){
  file.path(psc@out_dir, "pooled_bams", paste0(bam_group_name, "_meta.bam"))
}

#' Title
#'
#' @param psc
#' @param bams
#' @param name
#' @param pool_script
#' @param skip_peaksat
#'
#' @return
#' @export
#'
#' @examples
submit_meta_pool_jobs = function(psc, bam_groups, bam_group_names = names(bams), pool_script = get_pool_script()){
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
  cmds = sapply(names(bam_groups), function(nam){
    bams = bam_groups[[nam]]
    make_meta_pool_cmd(psc = psc, bams = bams, name = paste0(nam, "_meta"))
  })
  cmd_outs = lapply(cmds, system, intern = TRUE)
  if(psc@job_scheduler %in% c("SGE", "SLURM")){
    jids = sapply(cmd_outs, function(cmd_out){
      jid = capture_jid(cmd_out)
      jid
    })
  }else{
    jids = character()
  }
  PS_OPTIONS$PS_JOB_IDS = c(getOption("PS_JOB_IDS", character()), jids)
  names(jids) = get_meta_bam_file(psc, names(jids))
  jids
}
