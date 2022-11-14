library(peaksat)
library(magrittr)
library(seqsetvis)
library(data.table)

if(FALSE){
  bam_files = dir("bam_file_inputs/bam_files_atac/", pattern = "bam$", full.names = TRUE)
  bam_files = bam_files[!grepl("pool", bam_files)]
  
  groups = sapply(strsplit(basename(bam_files), "_"), function(x)x[1])
  grouped = split(bam_files, groups)
  
  nam = names(grouped)[1]
  for(nam in names(grouped)){
    to_pool = paste(grouped[[nam]], collapse = " ")
    pooled = file.path(dirname(grouped[[nam]][1]), paste0(nam, "_pooled.bam"))
    cmd = paste("  samtools merge", pooled, to_pool)
    cmd2 = paste("  samtools index", pooled)
    writeLines(c("#!/bin/bash", paste0("if [ ! -f ", pooled, " ]; then"), cmd, cmd2, "fi"), paste0(nam, "_merge.sh"))
    # system(cmd)
    # system(cmd2)
  }
  
  to_pool = paste(bam_files, collapse = " ")
  pooled = file.path(dirname(bam_files[1]), paste0("meta", "_pooled.bam"))
  cmd = paste("  samtools merge", pooled, to_pool)
  cmd2 = paste("  samtools index", pooled)
  writeLines(c("#!/bin/bash", paste0("if [ ! -f ", pooled, " ]; then"), cmd, cmd2, "fi"), paste0("meta", "_merge.sh"))
  # system(cmd)
  # system(cmd2)
}

# paste(pool_list$IgG$file, collapse = " ")
# paste(pool_list$Ikaros$file, collapse = " ")

bam_files = dir("bam_file_inputs/bam_files_atac/", pattern = "bam$", full.names = TRUE)

# pbmcapply::pbmclapply(bam_files[grepl("pool", bam_files)], function(f){
#   cap_reads(f, 200e6)
# })

# bam_files = bam_files[grepl("pool", bam_files)]

sp = strsplit(basename(bam_files), "_")

cells = sapply(sp, function(x)x[1])
marks = sapply(sp, function(x)x[2])
reps = sapply(sp, function(x)x[3])

table(cells)
table(marks)

# dir('bam_file_inputs/bam_files_cnr_PR_paper2')



chip_bams = bam_files

todo_df = data.frame(chip = chip_bams, stringsAsFactors = FALSE)
# todo_df$input = input_bam

stopifnot(file.exists(todo_df$chip))
# stopifnot(file.exists(todo_df$input))
stopifnot(file.exists(paste0(todo_df$chip, ".bai")))
# stopifnot(file.exists(paste0(todo_df$input, ".bai")))

todo_df = todo_df[grepl("200M", todo_df$chip), , drop = FALSE]

psc = peaksat_config(out_dir = "peaksat_outputs/peak_saturation_atac", stat = "qValue", is_PE = TRUE)
jobs = submit_peaksat_jobs(psc, todo_df$chip, NULL, await_completion = FALSE)
watch_jids()



cnt_dt = load_counts(psc)
plot_peak_saturation_lines(cnt_dt)

cnt_dt = load_counts(psc, min_signalValue = c(1, 5, 10))

cnt_dt$sample = old2new[cnt_dt$sample]

cnt_dt$signal_cutoff %>% table

cnt_dt[, c("cell", "treatment", "mark", "rep") := tstrsplit(sample, "_")]

debug(plot_peak_saturation_lines.facetted)
plot_peak_saturation_lines.facetted(cnt_dt[signal_cutoff == 5]) +
  coord_cartesian(xlim = c(0, 200e6))

class(cnt_dt$read_count)
class(cnt_dt$peak_count)

cnt_dt$read_count = as.numeric(cnt_dt$read_count)
cnt_dt$peak_count = as.numeric(cnt_dt$peak_count)

ko = cnt_dt[, .(read_count = max(read_count)), .(sample)][read_count > 205e6]$sample
cnt_dt = cnt_dt[!sample %in% ko]

ggplot(cnt_dt[signal_cutoff == 1], aes(x = read_count, y = peak_count, group = sample)) + 
  geom_path()

ggsave("tmp.png")

# show_available_results = function(psc){
#   res = dir(get_result_dir(psc))
#   res = res[res != "sub_logs"]
#   res = sub("^peak_saturation.", "", res)
#   res = sub(sub("results", "", basename(get_result_dir(psc))), "", res)
#   res
# }

show_available_results(psc)

old2new

a = names(old2new[1])
b = names(old2new[3])
min_signalValue = 1

options(mc.cores = 20)
conv_dt = load_convergence(psc, a, b, min_signalValue = 1)

conv_dt = load_convergence(psc, a, b, min_signalValue = 5)

conv_dt

plot_convergence_bars(conv_dt, old2new[a] , old2new[b])
ggsave("tmp1.pdf")
plot_convergence_trends(conv_dt, old2new[a] , old2new[b])
ggsave("tmp2.pdf")
