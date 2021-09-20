testthat::context("Basic")
# flipping viewGranges
library(peaksat)
library(testthat)
library(data.table)
library(ggplot2)

bam_files = dir(system.file("extdata", package = "ssvQC"), pattern = "^M.+.bam$", full.names = TRUE)
by_input = split(bam_files, ifelse(grepl("input", bam_files), "input", "chip"))

pc = peaksat_config(job_scheduler = "SGE", noModel = TRUE)

submit_peaksat_jobs(pc = pc,
                    treat_bams = by_input$chip,
                    ctrl_bams = by_input$input,
                    await_completion = FALSE)

cnt_dt = load_counts(pc)
plot_peak_saturation_lines(cnt_dt)
plot_peak_saturation_lines.facetted(cnt_dt)

meta_out = submit_meta_pool_jobs(pc = pc,
                                 bam_groups = by_input$chip,
                                 pool_script = bam_group_names = c("meta_CTCF"))

submit_peaksat_jobs(pc = pc,
                    treat_bams = names(meta_out),
                    ctrl_bams = by_input$input[1],
                    await_completion = FALSE,
                    hold_jid_map = meta_out)

