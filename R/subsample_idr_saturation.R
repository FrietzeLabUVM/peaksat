# library(data.table)
# library(magrittr)
# library(ggplot2)
# r1 = "data/CG_H4K5ac_MCF10-4.random1.bam"
# r2 = "data/CG_H4K5ac_MCF10-4.random2.bam"
# c1 = "/slipstream/home/conggao/ChIP_seq/MCF10_H4Kac/MCF10A_input/MCF10A_input_rep1.bam"
# c2 = c1
# #
# r1 = "data/CG_H4K8ac_MCF10-4.random1.bam"
# r2 = "data/CG_H4K8ac_MCF10-4.random2.bam"
# c1 = "/slipstream/home/conggao/ChIP_seq/MCF10_H4Kac/MCF10A_input/MCF10A_input_rep1.bam"
# c2 = c1
#
#
# stopifnot(file.exists(r1))
# stopifnot(file.exists(r2))
# stopifnot(file.exists(c1))
# stopifnot(file.exists(c2))
#
# # stat_arg = c("qVal_010" = "-q .01 -o results_qVal_010")
# #stat_arg = c("qVal_005" = "-q .005 -o results_qVal_005")
# # stat_arg = c("pVal_010" = "-p .01 -o results_pVal_010")
# stat_arg = c("pVal_050" = "-p .05 -o results_pVal_050")
# # stat_arg = c("pVal_100" = "-p .10 -o results_pVal_100")
#
# cmd = paste0(
#   "bash submit_subsample_peak.sh -t ",
#   r1,
#   " -c ",
#   c1,
#   " -n ",
#   paste(sep = "_", sub(".bam", "", basename(r1)), names(stat_arg), "vsInput"),
#   " ",
#   stat_arg,
#   ifelse(grepl("IgG", c1), " -pe ", ""))
# print(cmd)
# sub_out = system(cmd, intern = TRUE)
#
# sub_dt1 = data.table(out = sub_out, fraction = seq(.1, 1, .1))
# if(nrow(sub_dt1) > 0)
#   sub_dt1[, c("fraction", "i", "jid") := tstrsplit(out, " ", keep = c(1, 2, 5))]
#
# cmd = paste0(
#   "bash submit_subsample_peak.sh -t ",
#   r2,
#   " -c ",
#   c2,
#   " -n ",
#   paste(sep = "_", sub(".bam", "", basename(r2)), names(stat_arg), "vsInput"),
#   " ",
#   stat_arg,
#   ifelse(grepl("IgG", c2), " -pe ", ""))
# print(cmd)
# sub_out = system(cmd, intern = TRUE)
#
# sub_dt2 = data.table(out = sub_out)
# if(nrow(sub_dt2) > 0)
#   sub_dt2[, c("fraction", "i", "jid") := tstrsplit(out, " ", keep = c(1, 2, 5))]
#
# sub_dt1
# sub_dt2
#
# ###
# #plot peaks
#
# root_dir  = "~/R/peak_saturation/"
#
#
# load_counts = function(wd){
#   cn = c("name", "seed", "fraction", "stat", "cutoff", "PE", "input")
#
#   read_count_files = dir(wd, pattern = "read_count$", full.names = TRUE)
#   peak_count_files = dir(wd, pattern = "peak_count$", full.names = TRUE)
#
#
#   names(read_count_files) = sub(".bam.+", "", basename(read_count_files))
#   names(peak_count_files) = sub(".bam.+", "", basename(peak_count_files))
#
#   rc_dt = lapply(read_count_files, fread, col.names = c("read_count", cn)) %>%
#     rbindlist(idcol = "wd")
#   pc_dt = lapply(peak_count_files, fread, col.names = c("peak_count", cn)) %>%
#     rbindlist(idcol = "wd")
#
#
#   cnt_dt = merge(rc_dt, pc_dt[, .(name, peak_count)], by = "name")
#
#   cnt_dt[, c("fraction", "rep") := tstrsplit(name, "\\.", keep = 2:3)]
#   cnt_dt[, fraction := as.numeric(fraction)/100]
#
#   cnt_dt$sample = sub("peak_saturation.", "", basename(wd))
#
#   cnt_dt[order(read_count)]
# }
# # wds = dir("results_qVal_005", pattern = "^peak_saturation", full.names = TRUE)
# wds = dir(file.path(root_dir, "results_pVal_050"), pattern = "^peak_saturation.+H4K[58][Aa][Cc]", full.names = TRUE)
# wds = wds[!grepl("T47D", wds)]
# wds =wds[dir.exists(wds)]
#
# cnt_dt = rbindlist(lapply(wds, load_counts))
# cnt_dt[, c("source", "mark", "cell") := tstrsplit(sample, "[_\\.]", keep = 1:3)]
#
# p = ggplot(cnt_dt, aes(x = read_count, y = peak_count, color = sample)) +
#   # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
#   geom_point() +
#   stat_summary(fun = mean, geom = "line", aes(group = sample)) +
#   labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") +
#   facet_wrap(~mark, scales = "free") +
#   coord_cartesian(ylim = c(0, 500000))
# p
#
# #run idr
# peak_files = dir(dir("results_pVal_050/", pattern = "random", full.names = TRUE), pattern = "narrowPeak$", full.names = TRUE)
# peak_dt = data.table(file = peak_files)
# peak_dt[, seed := tstrsplit(basename(file), "[_\\.]", keep = 3)]
# peak_dt[, fraction := tstrsplit(basename(file), "\\.", keep = 2)]
# peak_dt[, group := tstrsplit(basename(dirname(file)), "[_\\.]", keep = c(6))]
# peak_dt[, mark := tstrsplit(basename(dirname(file)), "_", keep = c(3))]
# peak_dt[, .N, .(mark, group)]
#
# peak_lists = lapply(split(peak_dt, peak_dt$mark), function(x)split(x, x$group))
# p = 2
# i = 1
# j = 6
# pl = peak_lists[[1]]
# SCRIPT = "/slipstream/home/joeboyd/scripts/qc_framework/run_IDR.sh"
# for(i in seq_len(nrow(pl[[1]]))){
#   for(j in seq_len(nrow(pl[[2]]))){
#     f1 = pl[[1]]$file[i]
#     f2 = pl[[2]]$file[j]
#     vs_name = paste(pl[[1]]$group[i], pl[[2]]$group[j], sep = "_")
#     frac_name = paste(pl[[1]]$fraction[i], pl[[2]]$fraction[j], sep = "_")
#     mark_name = paste(pl[[1]]$mark[i], pl[[2]]$mark[j], sep = "_")
#     message(i, " ", j)
#     message(vs_name)
#     message(frac_name)
#     message(f1, " and ", f2)
#     dir.create("/slipstream/home/joeboyd/R/peak_saturation/results_pVal_050/idr", showWarnings = FALSE)
#
#     cmd = paste("qsub -cwd -o sub_logs -e sub_logs", SCRIPT, f1, f2, paste0("/slipstream/home/joeboyd/R/peak_saturation/results_pVal_050/idr/", mark_name, ".", vs_name, ".", frac_name))
#     system(cmd)
#   }
# }
#
# #plot idr
# idr_files = dir("results_pVal_050/idr/", pattern = "npeaks-aboveIDR", full.names = TRUE)
# fread(idr_files[1])
# idr_dt = data.table(file = idr_files)
# idr_dt[, mark := tstrsplit(basename(file), "_", keep = 1)]
# idr_dt[, frac_name := tstrsplit(basename(file), "\\.", keep = 3)]
# idr_dt[, fraction1 := tstrsplit(frac_name, "_", keep = 1)]
# idr_dt[, fraction2 := tstrsplit(frac_name, "_", keep = 2)]
# idr_dt[,fraction1 := as.factor(as.numeric(fraction1)/100)]
# idr_dt[,fraction2 := factor(as.numeric(fraction2)/100)]
#
# my_fread_idr_peaks = function(f, idr = .05){
#   dt = fread(f)
#   dt[V3 == idr]$V4
# }
#
# idr_dt[, idr_peaks := my_fread_idr_peaks(file), .(file)]
#
#
# library(ggplot2)
# ggplot(idr_dt, aes(x = fraction1, y = fraction2, fill = idr_peaks, label = idr_peaks)) +
#   geom_tile() +
#   facet_wrap(~mark)
#
# ggplot(idr_dt, aes(x = fraction1, y = fraction2, fill = idr_peaks, label = round(idr_peaks/1000, 1))) +
#   geom_text() +
#   facet_wrap(~mark)
#
# ggplot(idr_dt, aes(x = mark, y = idr_peaks)) +
#   geom_boxplot()
#
#
# # idr_dt[, frac_name := tstrsplit(basename(file), "\\.", keep = 3)]
