library(peaksat)

bam_files = dir("bam_file_inputs/bam_files_H4_input/", pattern = "bam$", full.names = TRUE)

bam_files.input = bam_files[grepl("input", basename(bam_files))]
bam_files.input = bam_files.input[grepl("AT1", bam_files.input)]
# bam_files.input = bam_files.input[10]

bam_files.chip = bam_files[!grepl("input", basename(bam_files))]

sp = strsplit(basename(bam_files), "_")

cells = sapply(sp, function(x)x[1])
marks = sapply(sp, function(x)x[2])
reps = sapply(sp, function(x)x[3])

table(cells)
table(marks)

for(in_f in bam_files.input){
  psc = peaksat_config(out_dir = paste0("peak_saturation_10A_H4ac_vs_", basename(in_f)))
  
  todo_df = data.frame(chip = bam_files.chip, stringsAsFactors = FALSE)
  todo_df$input = in_f
  
  jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)
}

watch_jids()

