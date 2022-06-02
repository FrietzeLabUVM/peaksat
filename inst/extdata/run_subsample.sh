#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=1                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=10000                        # Memory total in MB (for all cores)
#SBATCH -o sub_logs/peaksaturation_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e sub_logs/peaksaturation_%j.err                 # File to which STDERR will be written, including job ID


while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam) bam="$2"; shift ;;
        -f|--fraction) fraction="$2"; shift ;;
        -n|--name) name="$2"; shift ;;
        -s|--seed) seed="$2"; shift ;;
        -wd|--workdir) wd="$2"; shift ;;
        -st|--samtools) samtools_path="$2"; shift;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $samtools_path ]; then echo --samtools samtools path is required. exit 1; fi
if [ ! -f $samtools_path ]; then echo --samtools $samtools_path file not found. exit 1; fi


echo bam is $bam
echo fraction is $fraction
echo name is $name
echo seed is $seed
echo workdir is $wd


if [ -z $wd ]; then wd=$(pwd)/subsampled_bams; fi
mkdir -p $wd

prefix=$wd/${name}
samp_bam=${prefix}.bam
if [ -z $seed ]; then
  seed=$(date +%N)
fi

p_cmd="$stat"

if [ -z $input_bam ]; then
  c_cmd=""
else
  c_cmd="-c $input_bam"
fi

log=${prefix}.log

#if [ -f ${samp_bam}.read_count ] && [ -f ${samp_bam}.peak_count ]; then
#  echo Results are present, delete ${samp_bam}.read_count or ${samp_bam}.peak_count to allow rerun.
#  exit 0
#fi

if [ ! -f $samp_bam ]; then
  if [ $fraction = 1 ]; then
    ln -s $bam $samp_bam
    if [ -f ${bam}.bai ]; then
      ln -s ${bam}.bai ${samp_bam}.bai
    else
      $samtools_path index ${samp_bam}
    fi
    echo fraction is 1 >> $log
    echo "$samtools_path view -c $bam to ${samp_bam}.read_count" >> $log
    echo $($samtools_path view -c $bam) ${name} ${seed} ${fraction} ${stat} ${pe} ${input} > ${samp_bam}.read_count
  else
    echo "$samtools_path view -s ${seed}${fraction} -b $bam to $samp_bam" >> $log
    $samtools_path view -s ${seed}${fraction} -b $bam > $samp_bam
    echo "$samtools_path index $samp_bam" >> $log
    $samtools_path index $samp_bam
    echo "$samtools_path view -c $samp_bam to ${samp_bam}.read_count" >> $log
    echo $($samtools_path view -c $samp_bam) ${name} ${seed} ${fraction} ${stat} ${pe} ${input} > ${samp_bam}.read_count
  fi
fi
