# In the bash terminal run: 

# ========= USER VARIABLES =========
BASE_DIR="/home/yhezema/scratch/MAGIC"
ASSEMBLY_DIR="$BASE_DIR/Assembly/assemblies"
LOG_DIR="$BASE_DIR/Assembly/logs/rb_denovo"
REF_GENOME="$BASE_DIR/reference/Gmax_880_v6.0.fai"

# Create directories
mkdir -p "$ASSEMBLY_DIR/${parent}" "$LOG_DIR"
WORK="$ASSEMBLY_DIR/${parent}"
LOG="$LOG_DIR/ref_index.log"

# ========= MODULES =========
module load StdEnv/2023
module load samtools/1.20
module load bowtie2/2.5.4
module load htslib/1.19

# ========= STEP 0: Prepare Reference ========= (15 min)
echo "[Step 2] Decompressing and indexing reference started at [$(date)]" >> $LOG_DIR/ref_index.log

# 1. Decompress reference 
zcat $REF_GENOME > ${REF_GENOME%.gz} 2>> $LOG_DIR/ref_index.log 

# 2. Index for samtools
samtools faidx ${REF_GENOME%.gz} >> $LOG_DIR/ref_index.log 2>&1

# 3. Build Bowtie2 index 
bowtie2-build ${REF_GENOME%.gz} $WORK/ref_index >> $LOG_DIR/ref_index.log 2>&1

echo "[Step 2] Decompressing and indexing reference finished at [$(date)]" >> $LOG_DIR/ref_index.log
#=============================================================================
#-----------------------------------------------------------------------------
# ================================================================
# STEP 2: Map reads against reference & define blocks (9-14 hours)
# ================================================================
#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=map_%A_%a.log
#SBATCH --array=0-7

# List of parents
parents=(1355 Prosper RG11 RG22 RG23 RG46 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}


# To ensures script stops on any error
set -euo pipefail

# ========= USER VARIABLES =========
BASE_DIR="/home/yhezema/scratch/MAGIC"
TRIM_DIR="$BASE_DIR/Assembly/trimmed"
ASSEMBLY_DIR="$BASE_DIR/Assembly/assemblies"
LOG_DIR="$BASE_DIR/Assembly/logs/rb_denovo"
REF_GENOME_bt2="$ASSEMBLY_DIR/ref_bt2"
parent="RG46" # Change to the desired parent each time

# ========= Create directories =========
mkdir -p "$ASSEMBLY_DIR/${parent}" "$LOG_DIR"
WORK="$ASSEMBLY_DIR/${parent}"
LOG="$LOG_DIR/${parent}_map.log"

# ========= MODULES =========
module load StdEnv/2023
module load samtools/1.20
module load bedtools/2.31.0
module load bowtie2/2.5.4

# ========= define trimmed reads =========
R1_pe="${TRIM_DIR}/${parent}/${parent}_R1_trimmed_pe.fastq.gz"
R2_pe="${TRIM_DIR}/${parent}/${parent}_R2_trimmed_pe.fastq.gz"
R1_se="${TRIM_DIR}/${parent}/${parent}_R1_trimmed_se.fastq.gz"
R2_se="${TRIM_DIR}/${parent}/${parent}_R2_trimmed_se.fastq.gz"
SE_COMBINED="${R1_se},${R2_se}"

echo "Start Step 2 for $parent at [$(date)] " >> "$LOG"

# ========= STEP 2.1: Map Reads =========
echo "[Step 2] $parent mapping paired + unpaired reads with bowtie2 started at [$(date)]" >> "$LOG"

# Map reads
bowtie2 --fast-local -p 8 \
  -x "$REF_GENOME_bt2/ref_index" \
  -1 "$R1_pe" -2 "$R2_pe" \
  -U "$SE_COMBINED" \
  -S "$WORK/${parent}_aligned.sam" 2>> "$LOG"

echo "[Step 2.1] $parent mapping paired + unpaired reads with bowtie2 completed successfully at [$(date)]" >> "$LOG"
#=============================================================================
#-----------------------------------------------------------------------------
#=============================================================================
# Continue mapping
===================
#!/bin/bash
#SBATCH --job-name=sort_array
#SBATCH --time=3:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=sort_%j.log
#SBATCH --output=sort_%A_%a.log
#SBATCH --array=0-7

# To ensures the script stops on any error 
set -euo pipefail

# ========= USER VARIABLES =========
BASE_DIR="/home/yhezema/scratch/MAGIC"
ASSEMBLY_DIR="$BASE_DIR/Assembly/assemblies"
LOG_DIR="$BASE_DIR/Assembly/logs/rb_denovo"

# Array of parent IDs
parents=(1355 Prosper RG11 RG22 RG23 RG64 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# Create directories
mkdir -p "$ASSEMBLY_DIR/${parent}" "$LOG_DIR"
WORK="$ASSEMBLY_DIR/${parent}"
LOG="$LOG_DIR/${parent}_map.log"

# ========= MODULES =========
#module --force purge
module load StdEnv/2023
module load samtools/1.20

echo "Continue Step 2 for $parent at [$(date)] " >> "$LOG"

# ========= STEP 2.2: Sorting aligned Reads =========
# Converting sam to bam
echo "[Step 2.2] $parent BAM conversion started at [$(date)]" >> "$LOG"
samtools view -@8 -b "$WORK/${parent}_aligned.sam" > "$WORK/${parent}_aligned_unsorted.bam"

# Remove sam to save space
echo "Removing intermediate SAM file at [$(date)]" >> "$LOG"
rm "$WORK/${parent}_aligned.sam"

# Sort bam file
echo "[Step 2.2] $parent sorting started at [$(date)]" >> "$LOG"

samtools sort -@8 -m 12G \
  -T "$WORK/${parent}_tmp" \
  -O bam \
  "$WORK/${parent}_aligned_unsorted.bam" \
  > "$WORK/${parent}_mapped_all.bam" 2>> "$LOG"

echo "[Step 2.2] $parent sorting completed successfully at [$(date)]" >> "$LOG"

# Remove unsorted bam to save space

echo "Removing unsorted bam file at [$(date)]" >> "$LOG"

rm "$WORK/${parent}_aligned_unsorted.bam"

# ========= STEP 2.3: Index sorted aligned Reads =========
echo "[Step 2.3] $parent index started at [$(date)]" >> "$LOG"

samtools index "$WORK/${parent}_mapped_all.bam" 2>> "$LOG"

echo "[Step 2.3] $parent index completed successfully at [$(date)]" >> "$LOG"

# ========= STEP 2.3: Index sorted aligned Reads =========
echo "[Step 2.3] $parent index started at [$(date)]" >> "$LOG"

samtools index "$WORK/${parent}_mapped_all.bam" 2>> "$LOG"

echo "[Step 2.3] $parent index completed successfully at [$(date)]" >> "$LOG"

#=============================================================================









