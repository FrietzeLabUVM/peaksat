library(peaksat)
# run input saturation, inputs are matched
in_dir = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/"

run_merged = TRUE

if(run_merged){
  input_bams.all = "bam_file_inputs/bam_files_CTCF/MCF10A-AT1-CA1_input_merged.bam"
  chip_bams.all = "bam_file_inputs/bam_files_CTCF/MCF10A-AT1-CA1_CTCF_merged.bam"
}else{
  input_bams.all = dir(dir(in_dir, pattern = 'input', full.names = TRUE), pattern = "bam$", full.names = TRUE)  
  chip_bams.all = dir(dir(in_dir, pattern = 'CTCF', full.names = TRUE), pattern = "bam$", full.names = TRUE)
  chip_bams.all = chip_bams.all[!(grepl("pooled", chip_bams.all) | grepl("input", chip_bams.all))]
}

input_cells = sapply(strsplit(basename(input_bams.all), "_"), function(x)x[1])
input_bams.all = split(input_bams.all, input_cells)

chip_cells = sapply(strsplit(basename(chip_bams.all), "_"), function(x)x[1])
chip_bams.all = split(chip_bams.all, chip_cells)

bam = input_bams.all[[1]]
parse_Mreads = function(bam){
  sapply(strsplit(basename(bam), "\\."), function(sp){
    k = which(sp == "bam") - 1
    sp[k]
  })
}
parse_Mreads(input_bams.all[[1]])
parse_cells = function(bam){
  sapply(strsplit(basename(bam), "\\."), function(sp){
    # k = which(sp == "bam") - 1
    sp[1]
  })
}
parse_cells(input_bams.all[[1]])

make_name = function(bam){
  sub(".bam$", "", basename(bam))
}

stopifnot(setequal(names(chip_bams.all), names(input_bams.all)))
for(cl in names(chip_bams.all)){
  input_bams = input_bams.all[[cl]]
  chip_bams = chip_bams.all[[cl]]
  message("--", cl)
  for(input_bam in input_bams){
    # message(names(input_bam))
    
    nam = make_name(input_bam)
    message("----", nam)
    # message(paste(nam, collapse = "\n"))
    # }
    # todo_df = as.data.frame(expand.grid(chip_bams, input_bams))
    # colnames(todo_df) = c("chip", "input")
    # head(todo_df)
    
    todo_df = data.frame(chip = chip_bams)
    todo_df$input = input_bam
    # psc = peaksat_config(out_dir = paste0("peaksat_outputs.input_saturation_matched/peak_saturation_", nam, ".full_ChIP_input_saturation"))
    psc = peaksat_config(out_dir = paste0("peaksat_outputs.CTCF/peak_saturation_", nam), stat_value = .05)
    jobs = submit_peaksat_jobs(psc, todo_df$chip, todo_df$input, await_completion = FALSE)
  }
}
watch_jids()

seqsetvis::get_mapped_reads(input_bams)/1e6
seqsetvis::get_mapped_reads(chip_bams[grepl("pool", chip_bams)])/1e6
seqsetvis::get_mapped_reads(chip_bams[!grepl("pool", chip_bams)])/1e6
