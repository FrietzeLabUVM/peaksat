library(peaksat)
library(magrittr)
library(seqsetvis)
library(data.table)

if(FALSE){ #this created pooled bam files
  cfg_dt = fread("PR_config.csv")
  pool_list = split(cfg_dt, cfg_dt$mark)

  cmds = c(
    paste("samtools merge bam_file_inputs/bam_files_cnr_PR_paper2/WT_combo_H3K4me3_pooled.bam", paste(pool_list$H3K4me3$file, collapse = " ")),
    paste("samtools merge bam_file_inputs/bam_files_cnr_PR_paper2/WT_combo_IgG_pooled.bam", paste(pool_list$IgG$file, collapse = " ")),
    paste("samtools merge bam_file_inputs/bam_files_cnr_PR_paper2/WT_combo_Ikaros_pooled.bam", paste(pool_list$Ikaros$file, collapse = " "))
  )
  for(cmd in cmds){
    message(cmd)
    system(cmd)
  }

  system("samtools index bam_file_inputs/bam_files_cnr_PR_paper2/WT_combo_H3K4me3_pooled.bam")
  system("samtools index bam_file_inputs/bam_files_cnr_PR_paper2/WT_combo_IgG_pooled.bam")
  system("samtools index bam_file_inputs/bam_files_cnr_PR_paper2/WT_combo_Ikaros_pooled.bam")

  cfg_dt[, final_file := file.path("bam_file_inputs/bam_files_cnr_PR_paper2", paste0(paste(cell, treatment, mark, rep, sep = "_"), ".bam"))]
  cfg_dt[, cmd := paste("ln -s", file, final_file)]
  cfg_dt[, cmd.bai := paste("ln -s", paste0(file, ".bai"), paste0(final_file, ".bai"))]

  cfg_dt[, system(cmd), .(file)]
  cfg_dt[, system(cmd.bai), .(file)]
}

# paste(pool_list$IgG$file, collapse = " ")
# paste(pool_list$Ikaros$file, collapse = " ")

bam_files = dir("bam_file_inputs/bam_files_cnr_PR_paper2/", pattern = "bam$", full.names = TRUE)

sp = strsplit(basename(bam_files), "_")

cells = sapply(sp, function(x)x[1])
marks = sapply(sp, function(x)x[2])
reps = sapply(sp, function(x)x[3])

table(cells)
table(marks)

# dir('bam_file_inputs/bam_files_cnr_PR_paper2')

psc = peaksat_config(out_dir = "peaksat_outputs/peak_saturation_cnr2")

chip_bams = bam_files[!grepl("IgG", bam_files)]

old2new = basename(chip_bams)
old2new = sub(".bam$", "", old2new)
names(old2new) = sub(".bam$", "", basename(normalizePath(chip_bams)))

input_bam = bam_files[grepl("IgG_pool", bam_files)]

#match to input
todo_df = data.frame(chip = chip_bams, stringsAsFactors = FALSE)
todo_df$input = input_bam

stopifnot(file.exists(todo_df$chip))
stopifnot(file.exists(todo_df$input))
stopifnot(file.exists(paste0(todo_df$chip, ".bai")))
stopifnot(file.exists(paste0(todo_df$input, ".bai")))

#only certain bam files need to be run
k = grepl("frozen.+pool", todo_df$chip)
todo_df[k,]

jobs = submit_peaksat_jobs(psc, todo_df$chip[k], todo_df$input[k], await_completion = FALSE)
watch_jids()


cnt_dt = load_counts(psc, min_signalValue = c(1, 5, 10))
cnt_dt$sample = old2new[cnt_dt$sample]
cnt_dt$signal_cutoff %>% table
cnt_dt[, c("cell", "treatment", "mark", "rep") := tstrsplit(sample, "_")]

p_ik = plot_peak_saturation_lines.facetted(cnt_dt[signal_cutoff == 5 & mark == "Ikaros"]) +
  coord_cartesian(xlim = c(0, 40e6))
p_k4 = plot_peak_saturation_lines.facetted(cnt_dt[signal_cutoff == 5 & mark == "H3K4me3"]) +
  coord_cartesian(xlim = c(0, 40e6))
cowplot::plot_grid(p_ik, p_k4)
ggsave("plot_CNR.pdf", width = 11, height = 7)

plot_peak_saturation_lines.facetted(cnt_dt[signal_cutoff == 5]) +
  coord_cartesian(xlim = c(0, 40e6))
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
