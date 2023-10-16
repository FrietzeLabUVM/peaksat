#!/bin/bash
SCRIPTS=$(dirname "$(readlink -f "$0")")

stat_arg="-p .001"
pe_arg=""
js="SGE"
no_model=""
hold=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--treatment) chip_bam="$2"; shift ;;
        -c|--control) input_bam="$2"; shift ;;
        -n|--name) name="$2"; shift ;;
        -g|--genome) g_arg="$2"; shift ;;
        -m|--macs2) macs2_path="$2"; shift ;;
        -st|--samtools) samtools_path="$2"; shift;;
        -pe|--paired-end) pe_arg="-pe" ;;
        -p|--pval) stat_arg="-p $2"; shift;;
        -q|--qval) stat_arg="-q $2"; shift;;
        -wd|--workdir) wd="$2"; shift;;
        -o|--outpath) out="$2"; shift;;
        -js|--job-scheduler) js="$2"; shift;;
        -noModel|--noModel) no_model="--noModel";;
        -h|--hold_jid) hold_jid="$2"; shift;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done



if [ $js != "SGE" ] && [ $js != "SLURM" ] && [ $js != "bash" ]; then
  echo Job scheduler must be either SGE, SLURM, or bash.  Was ${js}.
  exit 1
fi

if [ -n "$hold_jid" ]; then
  if [ $js == "SGE" ]; then
    hold="-hold_jid $hold_jid"
  fi
  if [ $js == "SLURM" ]; then
    hold="-d afterok:$hold_jid"
  fi
fi

if [ -z $macs2_path ]; then echo --macs2 macs2 path is required. exit 1; fi
if [ ! -f $macs2_path ]; then echo --macs2 $macs2_path file not found. exit 1; fi
if [ -z $samtools_path ]; then echo --samtools samtools path is required. exit 1; fi
if [ ! -f $samtools_path ]; then echo --samtools $samtools_path file not found. exit 1; fi

if [ -z $g_arg ]; then
  g_arg="-g hs";
else
  g_arg="-g ${g_arg}";
fi
if [ -z $name ]; then name=$(basename $chip_bam .bam); fi
if [ -z $chip_bam ]; then echo need chip_bam as '-t|--treatment', quit!; exit 1; fi
if [ -z $input_bam ]; then
  c_cmd=""
else
  c_cmd="-c $input_bam"
fi


mkdir -p sub_logs

if [ -z "$wd" ]; then
  if [ -z $out ]; then
    out=$(pwd)
    wd=peak_saturation.${name}
  else
    wd=${out}/peak_saturation.${name}
  fi
fi

mkdir -p ${out}/sub_logs

for i in 2; do
  #for frac in 10 13 17 20 23 27 30 33 37 40 43 47 50 53 57 60 63 67 70 73 77 80 83 87 90 93 97; do
  for frac in 10 20 30 40 50 60 70 80 90; do
    if [ $js == SGE ]; then
      cmd="qsub $hold -cwd -o ${out}/sub_logs -e ${out}/sub_logs ${SCRIPTS}/run_subsample_peak.sh -m $macs2_path -st $samtools_path -t $chip_bam $c_cmd -f .${frac} -n subset.0${frac}.${i} -s ${i} -wd ${wd} $stat_arg $g_arg $pe_arg $no_model"
      echo $cmd
      if [ ! -f  ${wd}/subset.0${frac}.${i}_peaks.narrowPeak ]; then
        sub_out=$($cmd)
        echo $frac $i $sub_out
      fi
    elif [ $js == SLURM ]; then
      cmd="sbatch $hold ${SCRIPTS}/run_subsample_peak.sh -m $macs2_path -st $samtools_path -t $chip_bam $c_cmd -f .${frac} -n subset.0${frac}.${i} -s ${i} -wd ${wd} $stat_arg $g_arg $pe_arg $no_model"
      echo $cmd
      if [ ! -f  ${wd}/subset.0${frac}.${i}_peaks.narrowPeak ]; then
        sub_out=$($cmd)
        echo $frac $i $sub_out
      fi
    elif [ $js == bash ]; then
      bash ${SCRIPTS}/run_subsample_peak.sh -m $macs2_path -st $samtools_path -t $chip_bam $c_cmd -f .${frac} -n subset.0${frac}.${i} -s ${i} -wd ${wd} $stat_arg $g_arg $pe_arg $no_model
    else
      echo Unrecognized job_scheduler $js.
      exit 1
    fi

  done;
done

if [ ! -f  ${wd}/subset.100.1_peaks.narrowPeak ]; then
  if [ $js == SGE ]; then
    sub_out=$(qsub $hold -cwd -o ${out}/sub_logs -e ${out}/sub_logs ${SCRIPTS}/run_subsample_peak.sh -m $macs2_path -st $samtools_path -t $chip_bam $c_cmd -f 1 -n subset.100.1 -s 1 -wd ${wd} $stat_arg $g_arg $pe_arg $no_model)
    echo $frac $i $sub_out
  elif [ $js == SLURM ]; then
    sub_out=$(sbatch $hold ${SCRIPTS}/run_subsample_peak.sh -m $macs2_path -st $samtools_path -t $chip_bam $c_cmd -f 1 -n subset.100.1 -s 1 -wd ${wd} $stat_arg $g_arg $pe_arg $no_model)
    echo $frac $i $sub_out
  elif [ $js == bash ]; then
    bash ${SCRIPTS}/run_subsample_peak.sh -m $macs2_path -st $samtools_path -t $chip_bam $c_cmd -f 1 -n subset.100.1 -s 1 -wd ${wd} $stat_arg $g_arg $pe_arg $no_model
  else
    echo Unrecognized job_scheduler $js.
    exit 1
  fi
fi
