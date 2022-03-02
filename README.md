# peaksat
R package to do peak saturation analysis for ChIP-seq (and other) data.

## Setup Notes

Depends on macs2 and samtools being installed and executable.

So far this has only been tested on a Linux system but should work if macs2, samtools, and bash are present.

Will use either the Sun Grid Engine (SGE) or SLURM (This is a work in progress) job schedulers if they're present, otherwise will run all tasks in serial using standard bash and R's mclapply for parallelization.

## Installation

Installing from github is simple if you have devtools installed.

```
devtools::install_github("FrietzeLabUVM/peaksat")
```
