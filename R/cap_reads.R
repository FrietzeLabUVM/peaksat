get_mapped_reads = function(f){
  stats = Rsamtools::idxstatsBam(f)
  sum(stats[,3])
}

#' cap_reads
#'
#' create a subsampled version of bam_file that is approximately the number specified by max_reads
#'
#' @param bam_file an indexed bam file
#' @param max_reads target number of reads
#' @param ncores number of cores to use for samtools view operation
#'
#' @return path to subsampled file
#' @export
#'
#' @examples
#' cap_reads
cap_reads = function(bam_file, max_reads, ncores = 1, out_dir = dirname(bam_file)){
  # f = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_input_R1/MCF10A_input_R1.bam"
  if(!dir.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE)
  }
  res_file = function(f, max_reads){
    file.path(out_dir, sub(".bam", paste0(".", max_reads/1e6, "M.bam"), basename(f)))
  }
  out_f = res_file(bam_file, max_reads)
  message(bam_file, " to ", out_f)
  read_count = get_mapped_reads(bam_file)
  fraction = max_reads / read_count
  if(fraction >= 1){
    system(paste("ln -s", bam_file, out_f))
  }else{
    cmd=paste0("samtools view --threads ", ncores, " -s ", fraction, " -b ", bam_file, " > ", out_f)
    system(cmd)
    system(paste("samtools index", out_f))

  }
  out_f

}
