## Objective

# Demonstrate use of peaksat on CTCF.
# 
## Setup

# Load libraries and located files.

library(peaksat)
library(ggplot2)
library(cowplot)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(knitr)
options(mc.cores = 20)

out_dir = "paper_figures.110922.CTCF"

in_dir = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/"

input_bams.all = dir(dir(in_dir, pattern = 'input', full.names = TRUE), pattern = "bam$", full.names = TRUE)
input_cells = sapply(strsplit(basename(input_bams.all), "_"), function(x)x[1])
input_bams.all = split(input_bams.all, input_cells)

lapply(unlist(input_bams.all), seqsetvis::get_mapped_reads)

stat_value = .05
input_use_pval = FALSE
input_mode = "matched"

for(input_use_pval in c(FALSE)){
  for(input_mode in c("matched")){
    sub_dir = paste0("CTCF_", 
                     input_mode, "_", 
                     ifelse(input_use_pval, "pValue_", "qValue_"),
                     stat_value * 100)
    res_file = function(f){
      final_dir = file.path(out_dir, sub_dir)
      dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
      file.path(final_dir, f)
    }
    
    ## Input saturation
    
    #these paths are to outputs from running submit_peaksat_jobs
    # these scripts did the submitting
    # 
    # "run_peaksat.10a_CTCF.R" to 
    # peaksat_outputs.CTCF
    # 
    
    switch(input_mode,
           matched = {
             input_results = dir("peaksat_outputs.CTCF", full.names = TRUE)
             
             fig_height = 18
             input_mode_title = "vs cell line matched inputs"
           },
           merged = {
             stop()
             # input_results = dir("peaksat_outputs.input_saturation_metapool", pattern = "complete_ChIP", full.names = TRUE)
             # input_depths = sapply(strsplit(basename(input_results), "[\\._]"), function(x)x[5])
             # fig_height = 18
             # input_mode_title = "vs meta-pooled inputs"
           },
           stop("bad input_mode"))
    
    names(input_results) = sub("_input.+", "", sub("peak_saturation_", "", basename(input_results)))
    
    min_signalValue_todo = c(1, 2.5, 5)
    
    
    for(sig_cutoff in min_signalValue_todo){
      input_res = lapply(input_results, function(res_dir){
        message(res_dir)
        psc.in = peaksat_config(res_dir, stat = ifelse(input_use_pval, "pValue", "qValue"), stat_value = stat_value)
        load_counts(psc.in, min_signalValue = sig_cutoff)
      })
      dt.in = rbindlist(input_res, idcol = "input_depth")
      dt.in[, c("cell", "mark", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:3)]
      dt.in = dt.in[cell %in% c("MCF10A", "MCF10AT1", "MCF10CA1")]
      
      input_depth_lev = unique(dt.in$input_depth)
      dt.in$input_depth = factor(dt.in$input_depth, levels = input_depth_lev)
      
      dt.in = dt.in[order(read_count)]
      dt.in$read_depth = paste0(round(dt.in$read_count/1e6), "M")
      dt.in$read_depth = factor(dt.in$read_depth, levels = unique(dt.in$read_depth))
      
      dt.in[, peak_str := round(peak_count/1e3, 1)]
      dt.in$cell %>% table
      
      txt_size = 2.6
      
      dt.in$cell = factor(dt.in$cell, levels = c("MCF10A", "MCF10AT1", "MCF10CA1"))
      dt.in = dt.in[order(cell)]
      
      todo = expand.grid(
        rep = unique(dt.in$rep),
        cell = unique(dt.in$cell)
      )
      
      dt.in[, facet_x := paste(mark, rep)]
      dt.in$facet_x = factor(dt.in$facet_x, levels = unique(dt.in$facet_x))
      dt.in[, rep := sub("R", "rep", rep)]
      
      sel_mark = "CTCF"
      for(sel_mark in c("CTCF")){
        dt.in_sel = dt.in[mark == sel_mark]
        dt.in_sel = dt.in_sel[order(cell)]
        
        plot_title = paste("CTCF peak saturation")
        plot_subt = paste(sep = "\n",
                          paste("input", input_mode),
                          paste("signalValue >", sig_cutoff),
                          paste(ifelse(input_use_pval, "pValue", "qValue"), "<", stat_value)
        )
        p_sat = ggplot(dt.in_sel[signal_cutoff==sig_cutoff], 
                       aes(x = read_count, y = peak_count, group = sample, color = rep)) + 
          geom_path() +
          facet_grid(mark~cell, scales = "free_x") +
          scale_color_manual(values = c(rep1 = "#da0000ff", rep2 = "#0000faff")) +
          labs(x = "read count (M)", y = "peak count (k)", title = plot_title, subtitle = plot_subt) +
          scale_x_continuous(labels = function(x)x/1e6) +
          scale_y_continuous(labels = function(x)x/1e3) +
          cowplot::theme_cowplot() + 
          theme(axis.text = element_text(size = 6), legend.position = "bottom", title = element_text(size = 10), plot.subtitle = element_text(size = 8))
        
        plot(p_sat)
        ggsave(
          res_file(paste0("figCTCF_",
                          input_mode, "_",
                          "signalValue_",
                          sig_cutoff, "_", 
                          ifelse(input_use_pval, "pValue", "qValue"), ".lines.pdf")), 
          p_sat,
          width = 4.7, height = 3.15)
      }
    }
  }
}


for(input_use_pval in c(FALSE)){
  for(input_mode in c("matched")){
    sub_dir = paste0("CTCF_", 
                     input_mode, "_", 
                     ifelse(input_use_pval, "pValue_", "qValue_"),
                     stat_value * 100)
    res_file = function(f){
      final_dir = file.path(out_dir, sub_dir)
      dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
      file.path(final_dir, f)
    }
    
    ## Input saturation
    
    #these paths are to outputs from running submit_peaksat_jobs
    # these scripts did the submitting
    # 
    # "run_peaksat.10a_CTCF.R" to 
    # peaksat_outputs.CTCF
    # 
    
    switch(input_mode,
           matched = {
             input_results = dir("peaksat_outputs.CTCF", full.names = TRUE)
             
             fig_height = 18
             input_mode_title = "vs cell line matched inputs"
           },
           merged = {
             stop()
             # input_results = dir("peaksat_outputs.input_saturation_metapool", pattern = "complete_ChIP", full.names = TRUE)
             # input_depths = sapply(strsplit(basename(input_results), "[\\._]"), function(x)x[5])
             # fig_height = 18
             # input_mode_title = "vs meta-pooled inputs"
           },
           stop("bad input_mode"))
    
    names(input_results) = sub("_input.+", "", sub("peak_saturation_", "", basename(input_results)))
    
    min_signalValue_todo = c(1, 2.5, 5)
    
    
    for(sig_cutoff in min_signalValue_todo){
      input_res = lapply(input_results, function(res_dir){
        message(res_dir)
        psc.in = peaksat_config(res_dir, stat = ifelse(input_use_pval, "pValue", "qValue"), stat_value = stat_value)
        load_counts(psc.in, min_signalValue = sig_cutoff)
      })
      dt.in = rbindlist(input_res, idcol = "input_depth")
      dt.in[, c("cell", "mark", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:3)]
      dt.in = dt.in[!cell %in% c("MCF10A", "MCF10AT1", "MCF10CA1")]
      
      input_depth_lev = unique(dt.in$input_depth)
      dt.in$input_depth = factor(dt.in$input_depth, levels = input_depth_lev)
      
      dt.in = dt.in[order(read_count)]
      dt.in$read_depth = paste0(round(dt.in$read_count/1e6), "M")
      dt.in$read_depth = factor(dt.in$read_depth, levels = unique(dt.in$read_depth))
      
      dt.in[, peak_str := round(peak_count/1e3, 1)]
      dt.in$cell %>% table
      
      txt_size = 2.6
      
      # dt.in$cell = factor(dt.in$cell, levels = c("MCF10A", "MCF10AT1", "MCF10CA1"))
      dt.in = dt.in[order(cell)]
      
      todo = expand.grid(
        rep = unique(dt.in$rep),
        cell = unique(dt.in$cell)
      )
      
      dt.in[, facet_x := paste(mark, rep)]
      dt.in$facet_x = factor(dt.in$facet_x, levels = unique(dt.in$facet_x))
      dt.in[, rep := sub("R", "rep", rep)]
      
      sel_mark = "CTCF"
      for(sel_mark in c("CTCF")){
        dt.in_sel = dt.in[mark == sel_mark]
        dt.in_sel = dt.in_sel[order(cell)]
        
        plot_title = paste("CTCF peak saturation")
        plot_subt = paste(sep = "\n",
                          paste("input", input_mode),
                          paste("signalValue >", sig_cutoff),
                          paste(ifelse(input_use_pval, "pValue", "qValue"), "<", stat_value)
        )
        p_sat = ggplot(dt.in_sel[signal_cutoff==sig_cutoff], 
                       aes(x = read_count, y = peak_count, group = sample)) + 
          geom_path() +
          facet_grid(mark~cell, scales = "free_x") +
          # scale_color_manual(values = c(rep1 = "#da0000ff", rep2 = "#0000faff")) +
          labs(x = "read count (M)", y = "peak count (k)", title = plot_title, subtitle = plot_subt) +
          scale_x_continuous(labels = function(x)x/1e6) +
          scale_y_continuous(labels = function(x)x/1e3) +
          cowplot::theme_cowplot() + 
          theme(axis.text = element_text(size = 6), 
                legend.position = "bottom", 
                title = element_text(size = 8), 
                plot.subtitle = element_text(size = 8))
        
        plot(p_sat)
        ggsave(
          res_file(paste0("figCTCFmerge_",
                          input_mode, "_",
                          "signalValue_",
                          sig_cutoff, "_", 
                          ifelse(input_use_pval, "pValue", "qValue"), ".lines.pdf")), 
          p_sat,
          width = 2.5, height = 3.15)
      }
    }
  }
}
