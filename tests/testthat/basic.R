testthat::context("Basic")
# flipping viewGranges
library(peaksat)
library(testthat)
library(data.table)
library(ggplot2)

# bam_files = dir(system.file("extdata", package = "ssvQC"), pattern = "^M.+.bam$", full.names = TRUE)

#Key inputs are 2 parallel character vectors of bam files, 1 of ChIPs, 1 of inputs
#Parallel means that the Input for item 1 of ChIPs is item 1 of Input, 2 with 2, 3 with 3, etc.
#This script is meant to handle 1 mark at a time
#You could obviously wrap this in a loop to handle several marks simultaneously.

bam_files = dir(dir("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/", full.names = TRUE), pattern = "^M.+.bam$", full.names = TRUE)
by_input = split(bam_files, ifelse(grepl("input", bam_files), "input", "chip"))
by_input$chip = by_input$chip[1:6]
by_input$input = by_input$input[1:3]
by_input$input = rep(by_input$input, each = 2)

#pc = peaksat_config(job_scheduler = "SGE", noModel = TRUE, out_dir = "~/peaksat_dev")
pc = peaksat_config(job_scheduler = "SGE", out_dir = "~/peaksat_dev")

# Run peaksat on provided bams

ps1_jids = submit_peaksat_jobs(pc = pc,
                              treat_bams = by_input$chip,
                              ctrl_bams = by_input$input,
                              await_completion = FALSE)

# Merge all ChIPs into a single "meta" ChIP to get as many reads as possible.
meta_out = submit_meta_pool_jobs(pc = pc,
                                 bam_groups = by_input$chip,
                                 bam_group_names = c("combined_CTCF"))

# Run peaksat on meta ChIPs
ps2_jids = submit_peaksat_jobs(pc = pc,
                               treat_bams = names(meta_out),
                               ctrl_bams = by_input$input[1],
                               await_completion = FALSE,
                               hold_jid_map = meta_out)

watch_jids(c(ps1_jids, ps2_jids))

cnt_dt = load_counts(pc)
plot_peak_saturation_lines(cnt_dt)
plot_peak_saturation_lines.facetted(cnt_dt)
