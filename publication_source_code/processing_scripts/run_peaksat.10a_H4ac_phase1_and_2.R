library(peaksat)

bam_files = dir("bam_files_H4_ac_new_seq/", pattern = "bam$", full.names = TRUE)

sp = strsplit(basename(bam_files), "_")

cells = sapply(sp, function(x)x[1])
marks = sapply(sp, function(x)x[2])
reps = sapply(sp, function(x)x[3])

table(cells)
table(marks)



psc = peaksat_config(out_dir = "peak_saturation_10A_H4ac_seq1_and_seq2")

## run peaksat 
get_matching_input = function(bam_file){
  bam_files[marks == "input" & grepl("pooled", bam_files)]
}

chip_bams = c(
  bam_files[!grepl("input", bam_files)]
)

todo_df = data.frame(chip = chip_bams, stringsAsFactors = FALSE)
todo_df$input = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_input_R1/MCF10A_input_R1.bam"

jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)

watch_jids()
