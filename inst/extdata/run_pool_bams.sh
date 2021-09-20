#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=1                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=10000                        # Memory total in MB (for all cores)
#SBATCH -o sub_logs/pool_bams_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e sub_logs/pool_bams_%j.err                 # File to which STDERR will be written, including job ID

while [[ "$#" -gt 0 ]]; do
case $1 in
-b|--bams) bams="$2"; shift ;;
-n|--name) name="$2"; shift ;;
-st|--samtools) samtools_path="$2"; shift;;
-wd|--workdir) wd="$2"; shift ;;
*) echo "Unknown parameter passed: $1"; exit 1 ;;
esac
shift
done

if [ -z $samtools_path ]; then echo --samtools samtools path is required. exit 1; fi
if [ ! -f $samtools_path ]; then echo --samtools $samtools_path file not found. exit 1; fi

echo bams is $bams
echo name is $name
echo workdir is $wd
echo samtools_path is $samtools_path

if [ -z $wd ]; then wd=$(pwd)/pooled_bams; fi
mkdir -p $wd

out_bam=$wd/${name}.bam

log=$wd/${name}.log
echo Start $(date) > $log

if [ -f ${out_bam} ]; then
echo Results are present, delete ${out_bam} to allow rerun.
exit 0
fi

bams=${bams//,/" "}

merge_cmd="$samtools_path merge $out_bam $bams"
echo $merge_cmd >> $log
$merge_cmd
idx_cmd="$samtools_path index $out_bam"
echo $idx_cmd >> $log
$idx_cmd


echo Finish $(date) >> $log

