regBetween = function(x, reg_left, reg_capture, reg_right){
  reg_full = paste0("(?<=", reg_left, ")", reg_capture, "(?=", reg_right, ")")
  regmatches(x, regexpr(reg_full, x, perl = TRUE))
}

capture_jid = function(qsub_str){
  regBetween(qsub_str, "job ", "[0-9]+", " \\(")
}

watch_jids = function(hold_jids){
  pending = hold_jids
  last_len = length(pending)
  stat_str = ""
  while(length(pending) > 0){
    # message(length(pending))
    qstat_str = system("qstat", intern = TRUE)
    qstat_jid = regBetween(qstat_str, " ", "[0-9]+", " ")

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

qsub_and_wait = function(cmds){

  hold_jids = character()
  hold_jids = pbapply::pbsapply(cmds, function(cmd_sub){
    qsub_str = system(cmd_sub, intern = TRUE)
    jid = capture_jid(qsub_str)
    # message(paste0(jid, " "), appendLF = FALSE)
    # hold_jids = c(hold_jids, jid)
    jid
  })
  # for(cmd_sub in cmds){
  #
  # }
  # timestamp()
  watch_jids(hold_jids)
  # timestamp()
  hold_jids
}
