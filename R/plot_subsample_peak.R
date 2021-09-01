library(data.table)
library(ggplot2)
library(magrittr)
root_dir  = "~/R/peak_saturation_dev//"


load_counts = function(wd){
  cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
  
  read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
  peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)
  
  
  names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
  names(peak_count_files) = sub(".bam.+", "", basename(peak_count_files))
  
  rc_dt = lapply(read_count_files, fread, col.names = c("read_count", cn)) %>% 
    rbindlist(idcol = "wd")
  pc_dt = lapply(peak_count_files, fread, col.names = c("peak_count", cn)) %>% 
    rbindlist(idcol = "wd")
  
  
  cnt_dt = merge(rc_dt, pc_dt[, .(name, peak_count)], by = "name")
  
  cnt_dt[, c("fraction", "rep") := tstrsplit(name, "\\.", keep = 2:3)]
  cnt_dt[, fraction := as.numeric(fraction)/100]
  
  cnt_dt$sample = sub("peak_saturation.", "", basename(wd))
  
  cnt_dt[order(read_count)]
}
# wds = dir("results_qVal_005", pattern = "^peak_saturation", full.names = TRUE)
wds = dir(file.path(root_dir, "results_qVal_010"), pattern = "^peak_saturation.MCF7.+", full.names = TRUE)
wds = dir(file.path(root_dir, "results_qVal_010"), pattern = "^peak_saturation.+((H4K5[Aa][Cc])|(H4K8[Aa][Cc])).+", full.names = TRUE)
wds
wds = wds[!grepl("T47D", wds)]
wds =wds[dir.exists(wds)]

lens = sapply(wds, function(wd)length(dir(wd)))
# wds = wds[lens == max(lens)]

cnt_dt = rbindlist(lapply(wds, load_counts))

cnt_dt$sample %>% unique

cnt_dt$mark %>% table

cnt_dt[, mark := tstrsplit(sample, "_", keep = 2)]
cnt_dt[mark == "RUN", mark := tstrsplit(sample, "_", keep = 4)]
cnt_dt[grepl("h$", mark), mark := tstrsplit(sample, "_", keep = 3)]

cnt_dt$mark = toupper(cnt_dt$mark)

cnt_dt[, cell := tstrsplit(sample, "_", keep = 1)]
cnt_dt[cell == "RUN", cell := tstrsplit(sample, "_", keep = 3)]

cnt_dt$cell %>% table

cnt_dt$time = "na"
cnt_dt[, time := tstrsplit(sample, "_", keep = 2)]
cnt_dt[!grepl("h$", time), time := "na"]

# cnt_dt[, group := paste(cell, mark)]
# cnt_dt[, sample_split := sub("_qVal", "\nqVal", sample)]
# cnt_dt[, sample_split := sub("aligned_reads", "", sample_split)]
# cnt_dt[, sample_split := sub("CUT_RUN_", "", sample_split)]
# 
# cnt_dt[, source := "-"]
# cnt_dt[grepl("CUT_RUN", sample), source := "KQ_C&R"]
# cnt_dt[grepl("(^DF)|(^WT)", sample), source := "PR_C&R"]
# cnt_dt[grepl("(^MCF)|(^MDA)", sample), source := "ChIP"]
# table(cnt_dt$source)
# 
# # cnt_dt[, c("cell", "mark", "rep") := tstrsplit(sample_split, "[_\\.\n]", keep = 1:3)]
# 
# cnt_dt$rep %>% table
# cnt_dt[rep == "R1", rep := "rep1"]
# cnt_dt[rep == "R2", rep := "rep2"]

ggplot(cnt_dt, aes(x = read_count, y = peak_count, color = cell, group = sample)) +
  # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
  geom_point() +
  stat_summary(fun = mean, geom = "line", aes(group = sample)) +
  labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") +
  facet_wrap(~mark, scales = "free")

cnt_dt$sample = sub(".Align.+", "", cnt_dt$sample)
cnt_dt$sample = sub("_qVal_010+", "", cnt_dt$sample)
cnt_dt$sample = sub("_vsInput", "", cnt_dt$sample)
cnt_dt[, sample := gsub("_", "\n", sample)]

plots = lapply(unique(cnt_dt$mark), function(m){
  ggplot(cnt_dt[mark == m], aes(x = read_count, y = peak_count, color = cell, group = sample)) +
    # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
    geom_point() +
    stat_summary(fun = mean, geom = "line", aes(group = sample)) +
    labs(subtitle = "peaks called vs read count\nsubsampled in 10% increments", title = m, x= "reads (M)", y = "peaks (k)") +
    facet_wrap(~sample, scales = "free") +
    scale_x_continuous(labels = function(x)x/1e6) +
    scale_y_continuous(labels = function(x)x/1e3) +
    theme(panel.background = element_blank(), panel.grid = element_blank(), strip.text = element_text(size = 6), legend.position = "bottom")
})

fwrite(cnt_dt, "peak_saturation_H4K8AC_H4K5AC.csv")
pdf("peak_saturation_H4K8AC_H4K5AC.pdf", width = 8, height = 8.6)
plots[[1]]
plots[[2]]
dev.off()

p


p = ggplot(cnt_dt, aes(x = read_count, y = peak_count, color = source)) +
  # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
  geom_point() +
  stat_summary(fun = mean, geom = "line", aes(group = sample)) +
  labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") +
  facet_wrap(~sample_split, scales = "free")
p

cnt_dt[cell == "CG", cell := "4-cell-pool"]

p = ggplot(cnt_dt, aes(x = read_count, y = peak_count, color = cell)) +
  # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
  geom_point() +
  stat_summary(fun = mean, geom = "line", aes(group = sample)) +
  labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") +
  facet_wrap(~mark, scales = "free")
p


ggplot(cnt_dt[grepl("K16", sample)], aes(x = read_count, y = peak_count, color = rep)) +
  # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
  geom_point() +
  stat_summary(fun = mean, geom = "line", aes(group = sample)) +
  labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") #+
  #facet_wrap(~sample_split, scales = "free")


ggplot(cnt_dt[grepl("K20", sample)], aes(x = read_count, y = peak_count, color = rep)) +
  # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
  geom_point() +
  stat_summary(fun = mean, geom = "line", aes(group = sample)) +
  labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") #+
#facet_wrap(~sample_split, scales = "free")

sel_cnt_dt = cnt_dt[grepl("K4", sample)]

# ggplot(sel_cnt_dt[grepl("K4", sample)], aes(x = read_count, y = peak_count, color = cell)) +
#   # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
#   geom_point(data = sel_cnt_dt[rep == "pooled"], size = 2) +
#   geom_point(data = sel_cnt_dt[rep != "pooled"], size = 1) +
#   stat_summary(fun = mean, geom = "line", aes(group = sample)) +
#   labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") +
# facet_wrap(~source)

p

date_str = format(Sys.time(), "%y_%m_%d_%H%M")
ggsave(paste0("peak_saturation.", date_str, ".pdf"), p, width = 15, height = 9)
fwrite(cnt_dt, paste0("peak_saturation.", date_str, ".csv"))
