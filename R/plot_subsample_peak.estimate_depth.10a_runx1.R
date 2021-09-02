# library(data.table)
# library(ggplot2)
# library(magrittr)
# root_dir  = "~/R/peak_saturation_dev/"
#
# target_peaks = 2e4
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
#
# wds = c(
#   dir(file.path(root_dir, "results_qVal_010"), pattern = "^peak_saturation.AF.+", full.names = TRUE),
#   dir(file.path(root_dir, "results_qVal_010"), pattern = "^peak_saturation.10a_progression.+", full.names = TRUE)
# )
# wds
# wds =wds[dir.exists(wds)]
#
# lens = sapply(wds, function(wd)length(dir(wd)))
# # wds = wds[lens == max(lens)]
#
# cnt_dt = rbindlist(lapply(wds, load_counts))
#
# cnt_dt$sample %>% unique
#
# cnt_dt[, mark := tstrsplit(sample, "_", keep = 2)]
# cnt_dt$mark %>% table
#
# cnt_dt[mark == "progression", mark := "RUNX1"]
# cnt_dt[mark == toupper("PROGRESSION"), mark := "RUNX1"]
#
#
# cnt_dt$mark = toupper(cnt_dt$mark)
#
# cnt_dt[, cell := tstrsplit(sample, "_", keep = 1)]
# cnt_dt[cell == "10a", cell := "10a_combo"]
#
# cnt_dt$cell %>% table
#
# cnt_dt[, rep := tstrsplit(sample, "_", keep = 3)]
# cnt_dt$rep %>% table
#
# p = ggplot(cnt_dt, aes(x = read_count, y = peak_count, color = cell)) +
#   # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
#   geom_point() +
#   stat_summary(fun = mean, geom = "line", aes(group = sample)) +
#   labs(title = "Peaks called vs read count", subtitle = "subsampled in 10% increments") +
#   facet_wrap(~mark, scales = "free") +
#   scale_x_continuous(labels = function(x)paste0(x/1e6, "M")) +
#   scale_y_continuous(labels = function(x)paste0(x/1e3, "k"))
# p
#
# if(is.null(target_peaks)){
#   cnt_dt[, cutoff := .9*max(peak_count) , .(mark)]
# }else{
#   cnt_dt[, cutoff := target_peaks , .(mark)]
# }
#
# saturation_dt = cnt_dt[peak_count >= cutoff, .(saturation = min(read_count), read_count, peak_count), .(mark)]
# saturation_dt = saturation_dt[saturation == read_count]
# saturation_dt[, saturation := saturation / 1e6]
#
# saturation_dt[, slope := peak_count / read_count, .(mark)]
# saturation_dt[, read_count * slope]
#
# cnt_dt[, is_max := read_count == read_count[which.max(peak_count)], .(cell, mark, rep)]
# cnt_dt[, over_cutoff := peak_count > cutoff]
#
# min_dt = cnt_dt
#
# min_dt[, over_5m := read_count >= 5e6]
# min_dt = min_dt[over_5m == TRUE]
# min_dt[, is_min := read_count == min(read_count), .(cell, mark, rep)]
# min_dt = min_dt[is_min == TRUE]
# min_dt[, .(cell, mark, rep, read_count, peak_count)]
#
# max_dt = cnt_dt[is_max == TRUE | over_cutoff == TRUE]
# max_dt[, selected := read_count == min(read_count), .(cell, mark, rep) ]
# max_dt = max_dt[selected == TRUE]
#
# max_dt[, .(cell, mark, rep)]
#
#
#
# if(is.null(target_peaks)){
#   max_dt = merge(max_dt, saturation_dt[, .(mark, saturation_peak_count = peak_count)], by = "mark")
# }else{
#   max_dt = merge(max_dt, saturation_dt[, .(mark, saturation_peak_count = target_peaks)], by = "mark")
# }
#
#
# max_dt = merge(max_dt, min_dt[, .(zero_read_count = read_count, zero_peak_count = peak_count, cell, mark, rep)], by = c("cell", "mark", "rep"), all.x = TRUE)
#
# max_dt[is.na(zero_read_count), zero_read_count := 0]
# max_dt[is.na(zero_peak_count), zero_peak_count := 0]
#
#
# max_dt[, slope := (peak_count - zero_peak_count) / (read_count - zero_read_count)]
# max_dt[, yint := peak_count - (slope * read_count) ]
# max_dt[, saturation_read_count := (saturation_peak_count - yint) / slope]
#
# p_dt = max_dt[order(saturation_read_count), .(cell, mark, rep, saturation_read_count, peak_count, saturation_peak_count, yint)]
#
# ggplot(p_dt) +
#   geom_segment(aes(x = 0, y = yint, xend = saturation_read_count, yend = saturation_peak_count, color = rep)) +
#   facet_grid(cell~mark, scales = "free")
#
# p_estimate = ggplot(cnt_dt, aes(x = read_count, y = peak_count, color = rep)) +
#   # annotate("path", x= cnt_dt$read_count, y = cnt_dt$peak_count, color = rep("gray", nrow(cnt_dt)), group = cnt_dt$group) +
#   geom_point(data = cnt_dt[, .(peak_count == max(peak_count), peak_count, read_count), .(cell, mark, rep)][V1 == TRUE], alpha = 1) +
#   stat_summary(fun = mean, geom = "line", aes(group = sample), alpha = .5) +
#   geom_segment(data = p_dt, lty = 2,
#                aes(x = 0,
#                    y = yint,
#                    xend = saturation_read_count,
#                    yend = saturation_peak_count,
#                    color = rep)) +
#   labs(title = "Estimation of reads required for saturation", subtitle = "dashed line is extrapolation used for estimate\npoint is peak call with all reads") +
#   scale_x_continuous(labels = function(x)paste0(x/1e6, "M")) +
#   scale_y_continuous(labels = function(x)paste0(x/1e3, "k")) +
#   facet_wrap(cell~mark, scales = "free_x", ncol = 2)
# p_estimate
#
# p_dt[, M
#      _saturation_read_count := round(saturation_read_count / 1e6, 1)]
#
# out_dt = dcast(p_dt, cell+mark~rep, value.var = "M_saturation_read_count")
# out_dt = split(out_dt, out_dt$mark)
#
#
# date_str = format(Sys.time(), "%y_%m_%d_%H%M")
#
# openxlsx::write.xlsx(out_dt, paste0("peak_saturation_read_estimate.", date_str, ".xlsx"))
# ggsave(paste0("peak_saturation.", date_str, ".pdf"), p, width = 15, height = 9)
# ggsave(paste0("peak_saturation_read_estimate.", date_str, ".pdf"), p_estimate, width = 6, height = 8)
# fwrite(cnt_dt, paste0("peak_saturation.", date_str, ".csv"))
