library(peaksat)

bam_files.old = c(
  dir("bam_file_inputs/bam_files_H4_ac_old_seq", pattern = "bam$", full.names = TRUE),
  dir("bam_file_inputs/bam_files_H4_ac_new_seq/", pattern = "CA1a.+bam$", full.names = TRUE)
)
bam_files.old = bam_files.old[!grepl("_input_", bam_files.old)]

bam_files.new = dir("/slipstream/home/conggao/ChIP_seq/MCF10_H4Kac/", pattern = "bam$", full.names = TRUE)
bam_files.new = bam_files.new[grepl("_input_", bam_files.new)]

bam_files = c(bam_files.old, bam_files.new)

sp = strsplit(basename(bam_files), "[_\\.]")

cells = sapply(sp, function(x)x[1])
marks = sapply(sp, function(x)x[2])
reps = sapply(sp, function(x)x[3])

table(cells)
table(marks)
table(reps)

psc = peaksat_config(out_dir = "peaksat_outputs/peak_saturation_10A_H4ac_seq1_only.matched_inputs")

## run peaksat 
get_matching_input = function(bam_file){
  bam_files[marks == "input" & grepl("pooled", bam_files)]
}

chip_bams = c(
  bam_files[!grepl("_input_", bam_files)]
)
names(chip_bams) = cells[!grepl("_input_", bam_files)]

input_bams = c(
  bam_files[grepl("_input_", bam_files)]
)
names(input_bams) = cells[grepl("_input_", bam_files)]

todo_df = data.frame(chip = chip_bams, stringsAsFactors = FALSE)
todo_df$input = input_bams[names(chip_bams)]

todo_df

jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)

watch_jids()

seqsetvis::get_mapped_reads(input_bams)/1e6
seqsetvis::get_mapped_reads(chip_bams[grepl("pool", chip_bams)])/1e6
seqsetvis::get_mapped_reads(chip_bams[!grepl("pool", chip_bams)])/1e6

options(mc.cores = 10)
cnt_dt = load_counts(psc)

cnt_dt[, c("cell", "mark", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:3)]
cnt_dt[, .N, .(cell, mark, rep)]
