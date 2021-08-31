setClass("peaksat",
         representation = list(
           stat = "character",
           stat_value = "numeric",
           is_PE = "logical"
         ))

stopifnot(length(t0) == length(c0))

stat_arg = c("qVal_010" = "-q .01 -o results_qVal_010")
#stat_arg = c("qVal_005" = "-q .005 -o results_qVal_005")

for(i in seq_along(t0)){
  cmd = paste0(
    "bash submit_subsample_peak.sh -t ",
    t0[i],
    " -c ",
    c0[i],
    " -n ",
    paste(sep = "_", sub(".bam", "", basename(t0[i])), names(stat_arg), "vsInput"),
    " ",
    stat_arg,
    ifelse(grepl("IgG", c0[i]), " -pe ", ""))
  print(cmd)
  system(cmd)
}
