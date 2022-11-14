library(peaksat)
# run input saturation 
in_dir = "bam_file_inputs/bam_files_H4_input_new_seq"
bam_files = dir(in_dir, pattern = "bam$", full.names = TRUE)

## run peaksat 
# chip_bams = c(
#   bam_files[!grepl("input", basename(bam_files))]
# )
# chip_bams = chip_bams[grepl("25M", chip_bams)]
# names(chip_bams) = sub(".bam", "", sub("pooled.Aligned.sortedByCoord.out.", "", basename(chip_bams)))

chip_bams.all = dir("/slipstream/home/conggao/ChIP_seq/MCF10_H4Kac/", pattern = "bam$", full.names = TRUE)
chip_bams.all = chip_bams.all[!(grepl("pooled", chip_bams.all) | grepl("input", chip_bams.all))]
chip_cells = sapply(strsplit(basename(chip_bams.all), "_"), function(x)x[1])
chip_bams.all = split(chip_bams.all, chip_cells)
chip_bams = unlist(chip_bams.all)
names(chip_bams) = NULL

input_bams = c(
  bam_files[grepl("input", basename(bam_files))]
)
names(input_bams) = sub(".bam", "", sub("10A-AT1-CA1a-DCIS-", "", basename(input_bams)))
names(input_bams) = sub(".100M", "", names(input_bams))
names(input_bams)[!grepl("M$", names(input_bams))] = paste0(names(input_bams)[!grepl("M$", names(input_bams))], ".100M")

make_name = function(bam){
  names(bam) = sub(".bam", "", sub("10A-AT1-CA1a-DCIS-", "", basename(bam)))
  names(bam) = sub(".100M", "", names(bam))
  names(bam)[!grepl("M$", names(bam))] = paste0(names(bam)[!grepl("M$", names(bam))], ".100M")
  names(bam)
  
}

for(input_bam in input_bams){
  # message(names(input_bam))
  message(make_name(input_bam))
  nam = make_name(input_bam)
  # todo_df = as.data.frame(expand.grid(chip_bams, input_bams))
  # colnames(todo_df) = c("chip", "input")
  # head(todo_df)

  todo_df = data.frame(chip = chip_bams)
  todo_df$input = input_bam
  # psc = peaksat_config(out_dir = paste0("peaksat_outputs.input_saturation_metapool/peak_saturation_", nam, ".25M_ChIP_input_saturation"))
  psc = peaksat_config(out_dir = paste0("peaksat_outputs.input_saturation_metapool/peak_saturation_", nam, ".complete_ChIP_input_saturation"), stat = peaksat:::valid_stats$pValue, stat_value = .05)
  # jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)
  jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)
  
  psc2 = peaksat_config(out_dir = paste0("peaksat_outputs.input_saturation_metapool/peak_saturation_", nam, ".complete_ChIP_input_saturation"), stat = peaksat:::valid_stats$qValue, stat_value = .05)
  jobs2 = submit_peaksat_jobs(psc2, todo_df$chip, todo_df$input, await_completion = FALSE)
  
}

watch_jids()

seqsetvis::get_mapped_reads(input_bams)/1e6
seqsetvis::get_mapped_reads(chip_bams[grepl("pool", chip_bams)])/1e6
seqsetvis::get_mapped_reads(chip_bams[!grepl("pool", chip_bams)])/1e6
