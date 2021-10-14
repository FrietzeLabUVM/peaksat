get_submit_command = function(job_scheduler = valid_job_schedulers[1]){
  if(job_scheduler == "SGE"){
    "qsub -cwd"
  }else if(job_scheduler == "SLURM"){
    "sbatch"
  }else if(job_scheduler == "bash"){
    "bash"
  }else{
    stop("job_scheduler must be one of SGE, SLURM, or bash was ", job_scheduler)
  }
}

regBetween = function(x, reg_left, reg_capture, reg_right){
  reg_full = paste0("(?<=", reg_left, ")", reg_capture, "(?=", reg_right, ")")
  regmatches(x, regexpr(reg_full, x, perl = TRUE))
}

capture_jid = function(qsub_str, job_scheduler = "SGE"){
  if(job_scheduler == "SGE"){
    regBetween(qsub_str, "job ", "[0-9]+", " \\(")
  }else if(job_scheduler == "SLURM"){
    #qsub_str = "Submitted batch job 5701948"
    str = strsplit(qsub_str, " ")[[1]]
    str[length(str)]
  }else if(job_scheduler == "bash"){
    character()
  }else{
    stop("job_scheduler must be one of SGE or SLURM was ", job_scheduler)
  }
}

#' watch_jids
#'
#' @param hold_jids
#'
#' @return
#' @export
#'
#' @examples
watch_jids = function(hold_jids = getOption("PS_JOB_IDS", NULL), job_scheduler = "SGE"){
  if(is.null(hold_jids)){
    message("No hold_jids were given to watch.")
    return()
  }
  pending = hold_jids
  last_len = length(pending)
  stat_str = ""
  while(length(pending) > 0){
    # message(length(pending))

    if(job_scheduler == "SGE"){
      qstat_str = system("qstat", intern = TRUE)
      qstat_jid = regBetween(qstat_str, " ", "[0-9]+", " ")
    }else if(job_scheduler == "SLURM"){
      qstat_jid = system("squeue --format=%i", intern = TRUE)
    }else if(job_scheduler == "bash"){
      break
    }else{
      stop("job_scheduler must be one of SGE or SLURM was ", job_scheduler)
    }


    pending = intersect(pending, qstat_jid)
    cat(sprintf('\rpending jobs: %s ', paste0(length(pending), stat_str, paste(rep(" ", 60), collapse = ""))))
    if(length(pending) != last_len){
      stat_str = ""
      last_len = length(pending)
    }else{
      if(nchar(stat_str) > 25){
        stat_str = ""
      }else{
        stat_str = paste0(stat_str, "." )
      }
    }
    Sys.sleep(1)
  }
}

qsub_and_wait = function(cmds, job_scheduler = "SGE"){
  hold_jids = character()
  hold_jids = pbapply::pbsapply(cmds, function(cmd_sub){
    qsub_str = system(cmd_sub, intern = TRUE)
    jid = capture_jid(qsub_str, job_scheduler = job_scheduler)
    # message(paste0(jid, " "), appendLF = FALSE)
    # hold_jids = c(hold_jids, jid)
    jid
  })
  # for(cmd_sub in cmds){
  #
  # }
  # timestamp()
  watch_jids(hold_jids, job_scheduler = job_scheduler)
  # timestamp()
  hold_jids
}
