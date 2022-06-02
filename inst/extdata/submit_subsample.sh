#!/bin/bash
SCRIPTS=$(dirname "$(readlink -f "$0")")

js="SGE"
hold=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam) bam="$2"; shift ;;
        -n|--name) name="$2"; shift ;;
        -st|--samtools) samtools_path="$2"; shift;;
        -wd|--workdir) wd="$2"; shift;;
        -o|--outpath) out="$2"; shift;;
        -js|--job-scheduler) js="$2"; shift;;
        -h|--hold_jid) hold_jid="$2"; shift;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done



if [ $js != "SGE" ] && [ $js != "SLURM" ] && [ $js != "bash" ]; then
  echo Job scheduler must be either SGE, SLURM, or bash.  Was ${js}.
  exit 1
fi

if [ $js == "SGE" ]; then
  hold="-hold_jid $hold_jid"
fi
if [ $js == "SLURM" ]; then
  hold="-d afterok:$hold_jid"
fi

if [ -z $samtools_path ]; then samtools_path=$(which samtools); fi
if [ -z $samtools_path ]; then echo --samtools samtools path is required. exit 1; fi
if [ ! -f $samtools_path ]; then echo --samtools $samtools_path file not found. exit 1; fi

if [ -z $name ]; then name=$(basename $bam .bam); fi
if [ -z $bam ]; then echo need bam as '-b|--bam', quit!; exit 1; fi


mkdir -p sub_logs

if [ -z "$wd" ]; then
  if [ -z $out ]; then
    out=$(pwd)
    wd=downsampled_reads.${name}
  else
    wd=${out}/downsampled_reads.${name}
  fi
fi

mkdir -p ${out}/sub_logs


for i in 1; do
  #for frac in 10 13 17 20 23 27 30 33 37 40 43 47 50 53 57 60 63 67 70 73 77 80 83 87 90 93 97; do
  for frac in 10 20 30 40 50 60 70 80 90; do
    if [ $js == SGE ]; then
      cmd="qsub $hold -cwd -o ${out}/sub_logs -e ${out}/sub_logs ${SCRIPTS}/run_subsample.sh -st $samtools_path -b $bam -f .${frac} -n subset.0${frac}.${i} -s ${i} -wd ${wd}"
      if [ ! -f  ${wd}/subset.0${frac}.${i}_peaks.narrowPeak ]; then
        sub_out=$($cmd)
        echo $frac $i $sub_out
      fi
    elif [ $js == SLURM ]; then
      cmd="sbatch $hold ${SCRIPTS}/run_subsample.sh -st $samtools_path -b $bam -f .${frac} -n subset.0${frac}.${i} -s ${i} -wd ${wd}"
      if [ ! -f  ${wd}/subset.0${frac}.${i}_peaks.narrowPeak ]; then
        sub_out=$($cmd)
        echo $frac $i $sub_out
      fi
    elif [ $js == bash ]; then
      bash ${SCRIPTS}/run_subsample.sh -st $samtools_path -b $bam -f .${frac} -n subset.0${frac}.${i} -s ${i} -wd ${wd}
    else
      echo Unrecognized job_scheduler $js.
      exit 1
    fi

  done;
done


if [ $js == SGE ]; then
  sub_out=$(qsub $hold -cwd -o ${out}/sub_logs -e ${out}/sub_logs ${SCRIPTS}/run_subsample.sh -st $samtools_path -b $bam -f 1 -n subset.100.1 -s 1 -wd ${wd})
  echo $frac $i $sub_out
elif [ $js == SLURM ]; then
  sub_out=$(sbatch $hold ${SCRIPTS}/run_subsample.sh -st $samtools_path -b $bam -f 1 -n subset.100.1 -s 1 -wd ${wd})
  echo $frac $i $sub_out
elif [ $js == bash ]; then
  bash ${SCRIPTS}/run_subsample.sh -st $samtools_path -b $bam -f 1 -n subset.100.1 -s 1 -wd ${wd}
else
  echo Unrecognized job_scheduler $js.
  exit 1
fi

