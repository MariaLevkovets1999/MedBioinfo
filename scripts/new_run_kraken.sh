#!/bin/bash

#SBATCH --job-name=run_kraken            # Job name
##SBACH --output=/home/x_malev/applied_bioinfo/logs/%x_%A_%a.out  # Standard output
##SBATH --error=/home/x_malev/applied_bioinfo/logs/%x_%A_%a.err   # Standard error
#SBATCH --array=1-5                      # Array range (replace N with the actual number of pairs)
#SBATCH --cpus-per-task=2                # Number of CPU cores per task
#SBATCH --time=01:00:00                  # Time limit (hh:mm:ss)
#SBATCH --mem=90GB                       # Memory limit
#SBATCH --account=naiss2024-22-540       # Account name

# Set the data directory
DATA_DIR=/proj/applied_bioinformatics/users/x_malev/MedBioinfo/data/sra_fastq
OUT_DIR=/proj/applied_bioinformatics/users/x_malev/MedBioinfo/analysis/kraken
krona_dir=/proj/applied_bioinformatics/users/x_malev/MedBioinfo/analysis/krona
# Get the list of _1.fastq.gz files
file1=$(ls $DATA_DIR/*_1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
file2=${file1/_1.fastq.gz/_2.fastq.gz}
base_name=$(basename $file1 _1.fastq.gz)

# Run kraken2
srun --job-name=kracken2_${base_name} singularity exec -B /proj:/proj \
  /proj/applied_bioinformatics/common_data/kraken2.sif kraken2 \
  --db /proj/applied_bioinformatics/common_data/kraken_database/ \
  --threads 1 --paired --gzip-compressed \
  --report ${OUT_DIR}/${base_name}_report.txt \
  --output ${OUT_DIR}/${base_name}_output.txt \
  $file1 $file2

# Run bracken
srun --job-name=bracken_${base_name} singularity exec -B /proj/:/proj \
  /proj/applied_bioinformatics/common_data/kraken2.sif bracken \
  -d /proj/applied_bioinformatics/common_data/kraken_database/ \
  -i ${OUT_DIR}/${base_name}_bracken_report.txt \
  -o ${OUT_DIR}/${base_name}_bracken.txt

#needs to be -job_name kraken_run
sacct -P --format=JobID%15,JobName%18,ReqCPUS,ReqMem,Timelimit,State,ExitCode,Start,elapsedRAW,CPUTimeRAW,MaxRSS,NodeList -j ${SLURM_JOB_ID} | grep ERR > ${OUT_DIR}/kracken_bracken_sacct.txt

python /proj/applied_bioinformatics/tools/KrakenTools/kreport2krona.py -r ${OUT_DIR}/${base_name}_report_bracken_species.txt -o  ${krona_dir}/${base_name}.txt
#remove prefixes 
sed 's/__//g' ${krona_dir}/${base_name}.txt > ${krona_dir}/${base_name}.krona

srun --cpus-per-task=1 --time=00:00:10 --account=naiss2024-22-540 singularity exec -B /proj:/proj /proj/applied_bioinformatics/common_data/kraken2.sif ktImportText ${krona_dir}/${base_name}.krona -o ${krona_dir}/${base_name}.html

#srun --cpus-per-task=1 --time=00:00:10 --account=naiss2024-22-540 singularity exec -B /proj:/proj /proj/applied_bioinformatics/common_data/kraken2.sif ktImportText /proj/applied_bioinformatics/users/x_malev/MedBioinfo/analysis/krona/ERR6913113.krona -o /proj/applied_bioinformatics/users/x_malev/MedBioinfo/analysis/krona/ERR6913113.html 