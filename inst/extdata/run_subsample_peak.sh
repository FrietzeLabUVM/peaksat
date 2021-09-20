#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=1                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=10000                        # Memory total in MB (for all cores)
#SBATCH -o sub_logs/peaksaturation_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e sub_logs/peaksaturation_%j.err                 # File to which STDERR will be written, including job ID


se_cmd="--nomodel --extsize 180"
pe_cmd="--format BAMPE"
f_cmd=${se_cmd}
stat="-p .001"
pe=SE
input=noInput
no_model=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--treatment) chip_bam="$2"; shift ;;
        -c|--control) input_bam="$2"; input=vsInput; shift ;;
        -f|--fraction) fraction="$2"; shift ;;
        -n|--name) name="$2"; shift ;;
        -s|--seed) seed="$2"; shift ;;
        -wd|--workdir) wd="$2"; shift ;;
        -g|--genome) g="$2"; shift ;;
        -m|--macs2) macs2_path="$2"; shift ;;
        -st|--samtools) samtools_path="$2"; shift;;
        -pe|--paired-end) f_cmd=${pe_cmd}; pe=PE ;;
        -p|--pval) stat="-p $2"; shift;;
        -q|--qval) stat="-q $2"; shift;;
        -noModel|--noModel) no_model="--noModel";;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $macs2_path ]; then echo --macs2 macs2 path is required. exit 1; fi
if [ ! -f $macs2_path ]; then echo --macs2 $macs2_path file not found. exit 1; fi
if [ -z $samtools_path ]; then echo --samtools samtools path is required. exit 1; fi
if [ ! -f $samtools_path ]; then echo --samtools $samtools_path file not found. exit 1; fi


echo chip_bam is $chip_bam
echo input_bam is $input_bam
echo fraction is $fraction
echo name is $name
echo seed is $seed
echo workdir is $wd
echo genome is $g
echo pe is $pe
echo stat is $stat


if [ -z $g ]; then g=hs; fi
if [ -z $wd ]; then wd=$(pwd)/peak_saturation; fi
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

if [ -f ${samp_bam}.read_count ] && [ -f ${samp_bam}.peak_count ]; then
  echo Results are present, delete ${samp_bam}.read_count or ${samp_bam}.peak_count to allow rerun.
  exit 0
fi

if [ $fraction = 1 ]; then
  echo fraction is 1 >> $log
  echo "$samtools_path view -c $chip_bam to ${samp_bam}.read_count" >> $log
  echo $($samtools_path view -c $chip_bam) ${name} ${seed} ${fraction} ${stat} ${pe} ${input} > ${samp_bam}.read_count
  macs_cmd="$macs2_path callpeak -t $chip_bam ${c_cmd} -g $g -n $prefix ${p_cmd} ${f_cmd} ${no_model}"
  echo $macs_cmd >> $log
  $macs_cmd
else
  echo "$samtools_path view -s ${seed}${fraction} -b $chip_bam to $samp_bam" >> $log
  $samtools_path view -s ${seed}${fraction} -b $chip_bam > $samp_bam
  echo "$samtools_path index $samp_bam" >> $log
  $samtools_path index $samp_bam
  echo "$samtools_path view -c $samp_bam to ${samp_bam}.read_count" >> $log
  echo $($samtools_path view -c $samp_bam) ${name} ${seed} ${fraction} ${stat} ${pe} ${input} > ${samp_bam}.read_count
  macs_cmd="$macs2_path callpeak -t $samp_bam ${c_cmd} -g $g -n $prefix ${p_cmd} ${f_cmd} ${no_model}"
  echo $macs_cmd >> $log
  $macs_cmd
  if [ -f ${prefix}_peaks.narrowPeak ]; then
    rm ${samp_bam}
    rm ${samp_bam}.bai
  fi
fi

echo "cat ${prefix}_peaks.narrowPeak | wc -l to ${samp_bam}.peak_count" >> $log
echo $(cat ${prefix}_peaks.narrowPeak | wc -l)  ${name} ${seed} ${fraction} ${stat} ${pe} ${input} > ${samp_bam}.peak_count

