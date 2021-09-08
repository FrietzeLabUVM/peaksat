library(peaksat)
library(magrittr)
#locate files
bam_files = dir(dir("/slipstream/home/conggao/ChIP_seq/MCF10_Core_Targets", full.names = TRUE), full.names = TRUE, pattern = "\\.bam$")
bam_files = bam_files[file.exists(bam_files)]

dt = data.table(file = bam_files)
dt[, c("cell", "mark", "rep") := tstrsplit(basename(file), "[_\\.]", keep = 1:3)]
dt$cell %>% table

#match inputs
dt_inputs = dt[mark == "input"][, .(input_file = file, cell)]
file.exists(paste0(dt_inputs$input_file, ".bai"))
dt_inputs[, input_reads := ssvQC::get_mapped_reads(input_file), .(input_file)]
dt_inputs = dt_inputs[input_file != "/slipstream/home/conggao/ChIP_seq/MCF10_Core_Targets/MCF10A_dir/MCF10A_input_rep1.Aligned.sortedByCoord.out.bam"]

#the DCIS input is junk so replace with MCF10A
dt_inputs[cell == "MCF10DCIS"]$input_file = dt_inputs[cell == "MCF10A"]$input_file
dt_inputs

mark_grouped = split(dt, dt$mark)
mark_grouped = mark_grouped[setdiff(names(mark_grouped), "input")]

todo$H3K4me3
todo = lapply(mark_grouped, function(mark_dt){
  merge(mark_dt, dt_inputs, by = "cell")
})

pc = peaksat_config(out_dir = "~/peak_saturation/CG_MCF10_Core_Targets_2")
pc

lapply(todo, function(td){
  submit_peaksat_jobs(pc = pc, treat_bams = td$file, ctrl_bams = td$input_file, await_completion = FALSE)
})

cnt_dt = load_counts(pc)

#derive attributes and  cleanup sample (used for labels in plots)
cnt_dt[, c("cell", "mark", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:3)]
cnt_dt[, sample := sub(".Aligned.sortedByCoord.out_qValue_010", "", sample)]
cnt_dtl = split(cnt_dt, cnt_dt$mark)
names(cnt_dtl)

#you may want to a do a loop here but the data will need combed through manually I think.
m = "H3K27me3"
m = "H3K4me3"
#bad reps are poison (rep1 for DCIS K4me3 looks good actually)
plot_peak_saturation_lines(cnt_dtl[[m]])


m = "H3K27ac"
plot_peak_saturation_lines(cnt_dtl[[m]])
plot_peak_saturation_lines.facetted(cnt_dtl[[m]])

#linear regression does a decent overal job of estimating required read depth.
est_lin_res = estimate_depth.linear(cnt_dtl[[m]])
est_lin_res$estimates
est_lin_res$plots

#log regression does a better job on samples that have begun to saturate but a terrible job on those that haven't.
est_log_res = estimate_depth.log(cnt_dtl[[m]])
est_log_res$estimates
cowplot::plot_grid(plotlist = est_log_res$plots)

#ssvQC for fun
qdt = rbind(
  todo$H3K4me3[, .(file, cell, mark, rep)],
  dt_inputs[cell == "MCF10DCIS"][, .(file = input_file, cell = "MCF10A", mark = "input", rep = "rep1")]
)
qdt[, name := paste(cell, mark, rep, sep = "_")]
qdt[, name_split := gsub("_", "\n", name)]

np_files = dir(dir(get_result_dir(pc), pattern = "H3K4me3_rep", full.names = TRUE), pattern = "100.1_peaks.narrowPeak", full.names = TRUE)
np_dt = data.table(file = np_files)
np_dt[, name := tstrsplit(basename(dirname(file)), "\\.", keep = 2)]
np_dt[, name_split := gsub("_", "\n", name)]
np_dt[, c("cell", "mark", "rep") := tstrsplit(name, "_")]

library(ssvQC)
sqc = ssvQC(QcConfigFeatures(np_dt), QcConfigSignal(qdt))

sqc$signal_config$color_by = "name"
sqc$features_config$color_by = "name"
options(mc.cores = 20)
sqc = ssvQC.runAll(sqc)

sqc$plots$reads$All_signal
sqc$plots$features$count$All_features
sqc$plots$features$venn$All_features
sqc$plots$features$euler$All_features
sqc$plots$FRIP$reads_per_peak$All_features$All_signal
sqc$plots$signal$heatmaps$All_features$All_signal

sqc = ssvQC.selectClusters(sqc, cluster_numbers = c(5, 6))
sqc$plots$FRIP$per_peak$All_features$All_signal
sqc$plots$signal$heatmaps$All_features$All_signal
