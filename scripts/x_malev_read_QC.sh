#!/bin/bash

echo "1. Downloading the illumina sequencing data"
sqlite3 noheader -csv -batch /proj/applied_bioinformatics/common_data/sample_collab.db "select run_accession from sample_annot left join sample2bioinformatician s2b using(patient_code) where username='x_malev';" > /proj/applied_bioinformatics/users/x_malev/MedBioinfo/data/sra_fastq/x_malev_run_accessions.txt
cat /proj/applied_bioinformatics/users/x_malev/MedBioinfo/data/sra_fastq/x_malev_run_accessions.txt | srun --cpus-per-task=1 --time=00:30:00 --account=naiss2024-22-540  singularity exec /proj/applied_bioinformatics/common_data/meta.sif xargs fastq-dump --split-3 -gzip -O users/x_malev/MedBioinfo/data/sra_fastq/ --disable-multithreading --readids 

echo "2. Manipulate FASTQ files with seqkit"
singularity exec /proj/applied_bioinformatics/common_data/meta.sif seqkit -h
sacct --format=JobID,JobName%20,ReqCPUS,ReqMem,Timelimit,State,ExitCode,Start,elapsed,MaxRSS,NodeList,Account%15 
#Get general statistics: 
srun --cpus-per-task=20 --time=00:30:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit stats /mnt/data/ERR6913325.fastq
#Output: 
#file                        format  type  num_seqs      sum_len  min_len  avg_len  max_len
#/mnt/data/ERR6913325.fastq  FASTQ   DNA    912,111  217,026,247       70    237.9      302
#Now look for the number of unique reads 
srun --cpus-per-task=20 --time=00:30:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit rmdup -s /mnt/data/ERR6913325.fastq | singularity exec /proj/applied_bioinformatics/common_data/meta.sif seqkit stats
#Output: file  format  type  num_seqs      sum_len  min_len  avg_len  max_len
#-     FASTQ   DNA    844,345  201,464,992       70    238.6      302
echo "checking if sequencing kit adapters have been removed"
echo "testing on a random string like ATT"
srun --cpus-per-task=20 --time=00:30:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit grep -s -i -p ATT  /mnt/data/ERR6913325.fastq | grep -c "^@"
#result: 858530
echo "chekig with real kit adapters"
srun --cpus-per-task=20 --time=00:30:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit grep -s -i -p ATT  /mnt/data/ERR6913325.fastq | singularity exec /proj/applied_bioinformatics/common_data/meta.sif seqkit stats
#Output: 
#file  format  type  num_seqs      sum_len  min_len  avg_len  max_len
#-     FASTQ   DNA    858,530  206,724,613       70    240.8      302
echo 'Adaptor 1  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
srun --cpus-per-task=20 --time=00:30:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit grep -s -i -p  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  /mnt/data/ERR6913325.fastq | singularity exec /proj/applied_bioinformatics/common_data/meta.sif seqkit stats
echo 'Adaptor 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
srun --cpus-per-task=20 --time=00:30:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit grep -s -i -p AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  /mnt/data/ERR6913325.fastq | singularity exec /proj/applied_bioinformatics/common_data/meta.sif seqkit stats
echo "FASTQ files have already been trimmed of their sequencing kit adapters"

echo 'Quality control the raw sequencing FASTQ files with fastQC'
srun --cpus-per-task=2 --account=naiss2024-22-540 --time=00:10:00 singularity exec -B  /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif fastqc mnt/data/sra_fastq/ERR6913325.fastq  /mnt/data/sra_fastq/ERR6913296.fastq
srun --cpus-per-task=2 --account=naiss2024-22-540 --time=00:10:00 singularity exec -B  /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif xargs -a /mnt/data/sra_fastq/
x_malev_run_accessions.txt -I{} fastqc /mnt/data/sra_fastq/{}_1.fastq.gz  /mnt/data/sra_fastq/{}_2.fastq.gz
#merge files
srun --cpus-per-task=2 --account=naiss2024-22-540 --time=00:10:00 singularity exec -B  /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif xargs -a /mnt/data/sra_fastq/x_malev_run_accessions.txt -I{} flash -d /mnt/data/merged_pairs/{} /mnt/data/sra_fastq/{}_1.fastq.gz  /mnt/data/sra_fastq/{}_2.fastq.gz | tee -a /home/x_malev/applied_bioinfo/MedBioinfo/analysis/x_malev_flash.log
#or 
srun --cpus-per-task=2 --account=naiss2024-22-540 --time=00:10:00 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif xargs -a /mnt/data/sra_fastq/x_malev_run_accessions.txt -I{} sh -c '
echo 'Merging paired end reads'
flash /mnt/data/sra_fastq/{}_1.fastq.gz /mnt/data/sra_fastq/{}_2.fastq.gz \
&& mv out.extendedFrags.fastq /mnt/data/merged_pairs/{}_extendedFrags.fastq \
&& mv out.notCombined_1.fastq /mnt/data/merged_pairs/{}_notCombined_1.fastq \
&& mv out.notCombined_2.fastq /mnt/data/merged_pairs/{}_notCombined_2.fastq \
&& mv out.hist /mnt/data/merged_pairs/{}_.hist' | tee -a /home/x_malev/applied_bioinfo/MedBioinfo/analysis/x_malev_flash.log

echo 'made html files'

echo "Use read mapping to check for PhiX contamination"
echo "Download sequence data from the NCBI"
singularity exec /proj/applied_bioinformatics/common_data/meta.sif efetch -db nuccore -id NC_001422 -format fasta > /home/x_malev/applied_bioinfo/MedBioinfo/data/reference_seqs/PhiX_NC_001422.fna
srun --cpus-per-task=1 --time=00:30:00 --account=naiss2024-22-540 singularity exec /proj/applied_bioinformatics/common_data/meta.sif bowtie2-build -f ./reference_seqs/PhiX_NC_001422.fna ./bowtie2_DBs/PhiX_bowtie2_DB
srun  --cpus-per-task=8 --time=03:00:00 --account=naiss2024-22-540 singularity exec /proj/applied_bioinformatics/common_data/meta.sif bowtie2 -x ./data/bowtie2_DBs/PhiX_bowtie2_DB -U ./data/merged_pairs/ERR*extendedFrags.fastq -S ./analysis/bowtie/x_malev_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_malev_bowtie2_PhiX.log

#Day 1 exercise 
srun --cpus-per-task=2 --time=00:30:00 --account=naiss2024-22-
540 \bash -c "zcat /proj/applied_bioinformatics/common_data/refseq_viral_split/*.gz | \/proj/applied
_bioinformatics/tools/ncbi-blast-2.15.0+-src/makeblastdb \-out /home/x_malev/applied_bioinfo/MedBioi
nfo/data/blast_db  \-title 'RefSeq viral genomes' \-dbtype nucl"

srun --cpus-per-task=2 --time=00:05:00 --account=naiss2024-22-540 singularity exec -B /home/x_malev/applied_bioinfo/MedBioinfo/data:/mnt/data /proj/applied_bioinformatics/common_data/meta.sif seqkit sample -n 1000 /mnt/dat
a/sra_fastq/ERR6913114.fastq > ERR6913114_subet_1000.fastq

#Kraken
singularity exec /proj/applied_bioinformatics/common_data/kraken2.sif kraken2 -h

#pie chart 
