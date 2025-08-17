# ===============================
# STEP 1: Trimming with fastp
# ===============================
#!/bin/bash
#SBATCH --job-name=fastp_trim_%A_%a
#SBATCH --time=3:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7

# List of parents
parents=(1355 Prosper RG11 RG22 RG23 RG46 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# Base directory
BASE_DIR="/project/def-meskanda/yhezema/MAGIC"

# Create output directories
mkdir -p $BASE_DIR/Assembly/trimmed/${parent}
mkdir -p $BASE_DIR/Assembly/logs/fastp_reports

# Load fastp
module load StdEnv/2023 fastp/0.24.0

echo "[$(date)] Starting fastp trimming for $parent" >> $BASE_DIR/Assembly/logs/fastp_trimming_log.txt

# Run fastp
fastp \
  -i $BASE_DIR/parents/${parent}/NS.1509.002.D702---D502.${parent}_R1.fastq.gz \
  -I $BASE_DIR/parents/${parent}/NS.1509.002.D702---D502.${parent}_R2.fastq.gz \
  -o $BASE_DIR/Assembly/trimmed/${parent}/${parent}_R1_trimmed_pe.fastq.gz \
  -O $BASE_DIR/Assembly/trimmed/${parent}/${parent}_R2_trimmed_pe.fastq.gz \
  --unpaired1 $BASE_DIR/Assembly/trimmed/${parent}/${parent}_R1_trimmed_se.fastq.gz \
  --unpaired2 $BASE_DIR/Assembly/trimmed/${parent}/${parent}_R2_trimmed_se.fastq.gz \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --trim_poly_x \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 50 \
  --thread 4 \
  --html $BASE_DIR/Assembly/qc/fastp_reports/${parent}_report.html \
  --json $BASE_DIR/Assembly/qc/fastp_reports/${parent}_report.json \
  >> $BASE_DIR/Assembly/logs/fastp_trimming_log.txt 2>> $BASE_DIR/Assembly/logs/fastp_errors.log

echo "fastp trimming for $parent finished at $(date)" >> $BASE_DIR/Assembly/logs/fastp_trimming_log.txt
#-------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# QC with Fastqc   
#------------------------
#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7

# Define parent list
PARENTS=(1355 Prosper RG11 RG22 RG23 RG46 X790 07-78)
PARENT=${PARENTS[$SLURM_ARRAY_TASK_ID]}

# Load fastqc
module load fastqc/0.12.1 StdEnv/2023


echo "Starting processing ${PARENT} QC_Trimmed at $(date)" >> MAGIC/Assembly/logs/qc_trimmed_log.txt

# Run fastqc
fastqc \
  MAGIC/Assembly/trimmed/${PARENT}/${PARENT}_R1_trimmed_pe.fastq.gz \
  MAGIC/Assembly/trimmed/${PARENT}/${PARENT}_R1_trimmed_se.fastq.gz \
  MAGIC/Assembly/trimmed/${PARENT}/${PARENT}_R2_trimmed_pe.fastq.gz \
  MAGIC/Assembly/trimmed/${PARENT}/${PARENT}_R2_trimmed_se.fastq.gz \
  -o MAGIC/Assembly/qc/trimmed/ \
  2>> MAGIC/Assembly/logs/qc_trimmed_errors_log.txt

echo "Completed ${PARENT} QC at $(date)" >> MAGIC/Assembly/logs/qc_trimmed_log.txt
#=============================================================================
