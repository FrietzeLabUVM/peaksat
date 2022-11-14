library(peaksat)

bam_files = dir("/slipstream/home/conggao/ChIP_seq/MCF10_H4Kac/", pattern = "bam$", full.names = TRUE)

sp = strsplit(basename(bam_files), "[_\\.]")

cells = sapply(sp, function(x)x[1])
marks = sapply(sp, function(x)x[2])
reps = sapply(sp, function(x)x[3])

table(cells)
table(marks)
table(reps)

psc = peaksat_config(out_dir = "peaksat_outputs/peak_saturation_10A_H4ac_seq1_2_and_3.matched_inputs")

## run peaksat 
get_matching_input = function(bam_file){
  bam_files[marks == "input" & grepl("pooled", bam_files)]
}

chip_bams = c(
  bam_files[!grepl("input", bam_files)]
)
names(chip_bams) = cells[!grepl("input", bam_files)]

input_bams = c(
  bam_files[grepl("input", bam_files)]
)
names(input_bams) = cells[grepl("input", bam_files)]

todo_df = data.frame(chip = chip_bams, stringsAsFactors = FALSE)
todo_df$input = input_bams[names(chip_bams)]

jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)

watch_jids()

seqsetvis::get_mapped_reads(input_bams)/1e6
seqsetvis::get_mapped_reads(chip_bams[grepl("pool", chip_bams)])/1e6
seqsetvis::get_mapped_reads(chip_bams[!grepl("pool", chip_bams)])/1e6
