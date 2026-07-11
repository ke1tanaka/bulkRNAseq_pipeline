#!/bin/bash

## Maximum job time in Days-Hours:Minutes:Seconds
#SBATCH --time=1-00:00:00
## CPUs requested for each "task"; in simplest case the total number of used CPUs
#SBATCH --cpus-per-task=8
## Total memory; can also be expressed as --mem-per-cpu
#SBATCH --mem=120g

# load module
module reset
module load FastQC/0.12.1-Java-11

## QC check
fastqc -t 8 *.fastq.gz

## Load module
module reset
module load MultiQC/1.10.1-foss-2020b-Python-3.8.6

## MultiQC
multiqc -n multiqc_report_raw_reads.html .

# load module
module reset
module load fastp/0.23.2-GCCcore-10.2.0

## Fastp
## Loop through all *_1.fastq.gz files in the current directory
for file in *1_001.fastq.gz; do
    # Define the base name by removing the _1.fastq.gz suffix
    base=${file%1_001.fastq.gz}
    
    # Construct the corresponding pair file name
    pair="${base}2_001.fastq.gz"
    
    # Check if the pair file exists
    if [[ -f "$pair" ]]; then
        echo "Processing $file and $pair"
        
        # Run fastp on the paired files
        fastp -i "$file" -I "$pair" -o "${base}1_val_1.fq.gz" -O "${base}2_val_2.fq.gz" -l 50 --detect_adapter_for_pe -D --html "${base}_fastp_report.html" --thread 8
    else
        echo "Pair file for $file does not exist, skipping..."
    fi
done

# load module
module reset
module load FastQC/0.12.1-Java-11

## QC check
fastqc -t 8 *1_val_1.fq.gz
fastqc -t 8 *2_val_2.fq.gz

## Load module
module reset
module load MultiQC/1.10.1-foss-2020b-Python-3.8.6

## MultiQC
multiqc -n multiqc_report_trimmed_reads.html .

## Load module
module reset
module load STAR/2.7.11a-GCC-12.2.0

## Making STAR index
if [[ ! -f "STAR_index/SAindex" ]]; then
  mkdir STAR_index
  STAR \
  --runMode genomeGenerate \
  --genomeDir STAR_index \
  --runThreadN 8 \
  --genomeFastaFiles ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
fi

## Mapping by STAR
if [[ ! -f "STAR_output" ]]; then
  mkdir STAR_output
fi

for file in *1_val_1.fq.gz; do
    # Define the base name by removing the _1.fastq.gz suffix
    base=${file%1_val_1.fq.gz}
    
    # Construct the corresponding pair file name
    pair="${base}2_val_2.fq.gz"
    
    # Running STAR
    STAR \
      --genomeDir STAR_index \
      --runThreadN 8 \
      --outFileNamePrefix STAR_output/${base}_ \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn $file $pair \
        --readFilesCommand gunzip -c
done


## Remove duplicate reads using Picard
# load module
module reset
module load picard/2.25.6-Java-11

for file in STAR_output/*.bam; do
    # Define the base name by removing the .bam extension
    base=${file%.bam}
    
    # Running picard to remove duplicate reads
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      REMOVE_DUPLICATES=true \
      I=${file} \
      O=${base}_dedup.bam \
      M=${base}_dup_metrics.txt
done

## Load module
module reset
module load SAMtools/1.21-GCC-12.2.0
module load deepTools/3.5.5-foss-2022b

## Running deeptools to generate bigwig files
## Define the input and output directories
input_dir="STAR_output"
output_dir="deeptools_output"

## Index bam file
for bam_file in "$input_dir"/*_dedup.bam; do
   samtools index $bam_file
done

## Create the output directory if it doesn't exist
mkdir -p "$output_dir"

## Loop through all BAM files in the input directory
for bam_file in "$input_dir"/*_dedup.bam; do
  ## Extract the base name of the file (without path and extension)
  base_name=$(basename "$bam_file" .bam)

  ## Define the output file path
  output_file="$output_dir/${base_name}.bw"

  ## Run bamCoverage
  bamCoverage -b "$bam_file" -o "$output_file" --outFileFormat bigwig --numberOfProcessors 8
done

## Load module
module reset
module load Subread/2.0.3-GCC-10.2.0

## Run featureCounts
featureCounts -T 8 -p --countReadPairs -a ref/gencode.v49.primary_assembly.annotation.gtf -o count.out STAR_output/*_dedup.bam
