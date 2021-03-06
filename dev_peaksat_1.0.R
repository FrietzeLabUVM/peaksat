library(peaksat)
library(data.table)
library(ggplot2)

#Key inputs are 2 parallel character vectors of bam files, 1 of ChIPs, 1 of inputs
#Parallel means that the Input for item 1 of ChIPs is item 1 of Input, 2 with 2, 3 with 3, etc.
#This script is meant to handle 1 mark at a time
#You could obviously wrap this in a loop to handle several marks simultaneously.

bam_files = dir(dir("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/", full.names = TRUE), pattern = "^M.+.bam$", full.names = TRUE)
by_input = split(bam_files, ifelse(grepl("input", bam_files), "input", "chip"))
by_input$chip = by_input$chip[1:6]
by_input$input = by_input$input[1:3]
by_input$input = rep(by_input$input, each = 2)

#verify chip matches input
stopifnot(length(by_input$chip) == length(by_input$input))
check_df = data.frame(chip = by_input$chip, input = by_input$input)
check_df

psc = peaksat_config(job_scheduler = "SGE", out_dir = "~/peaksat_dev")

# Run peaksat on provided bams
ps1_jids = submit_peaksat_jobs(psc = psc,
                               treat_bams = by_input$chip,
                               ctrl_bams = by_input$input,
                               await_completion = FALSE)

# Merge all ChIPs into a single "meta" ChIP to get as many reads as possible.
meta_out = submit_meta_pool_jobs(psc = psc,
                                 bam_groups = by_input$chip,
                                 bam_group_names = c("combined_CTCF"))

# Run peaksat on meta ChIPs
ps2_jids = submit_peaksat_jobs(psc = psc,
                               treat_bams = names(meta_out),
                               ctrl_bams = by_input$input[1],
                               await_completion = FALSE,
                               hold_jid_map = meta_out)

#Create per cell pools
cells = sapply(strsplit(basename(by_input$chip), "_"), function(x)x[1])
by_cell_chip = split(by_input$chip, cells)
names(by_cell_chip) = paste0(names(by_cell_chip), "_pooled")
pooled_out = submit_meta_pool_jobs(psc = psc,
                                   bam_groups = by_cell_chip,
                                   bam_group_names = names(by_cell_chip))

by_cell_input = unique(by_input$input)

data.frame(chip = names(pooled_out), input = by_cell_input)

# Run peaksat on pooled ChIPs
ps3_jids = submit_peaksat_jobs(psc = psc,
                               treat_bams = names(pooled_out),
                               ctrl_bams = by_cell_input,
                               await_completion = FALSE,
                               hold_jid_map = pooled_out)



# It is critical that you wait for jobs to finish before loading the results
# You can capture the jids from each submit function or use value of PS_OPTIONS$PS_JOB_IDS which stores all jobs submitted this session
# watch_jids(c(ps1_jids, ps2_jids, ps3_jids))
watch_jids(PS_OPTIONS$PS_JOB_IDS)

cnt_dt = load_counts(psc)
plot_peak_saturation_lines(cnt_dt)
plot_peak_saturation_lines.facetted(cnt_dt)
