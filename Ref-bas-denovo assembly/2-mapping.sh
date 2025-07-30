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
#=============================================================================

# ================================================================
# STEP 2: Map reads against reference & define blocks (9-14 hours)
# ================================================================
#!/bin/bash
#SBATCH --job-name=map_RG11 # Change the name to the desired parent each time
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=map_%j.log

# To ensures script stops on any error
set -euo pipefail

# ========= USER VARIABLES =========
BASE_DIR="/home/yhezema/scratch/MAGIC"
TRIM_DIR="$BASE_DIR/Assembly/trimmed"
ASSEMBLY_DIR="$BASE_DIR/Assembly/assemblies"
LOG_DIR="$BASE_DIR/Assembly/logs/rb_denovo"
REF_GENOME_bt2="$ASSEMBLY_DIR/ref_bt2"
parent="RG11" # Change to the desired parent each time

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
# Continue mapping (30 min - 1 hour)
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

#=============================================================================
#-----------------------------------------------------------------------------
#=============================================================================
# Defining blocks and superblocks - all parents (4-5 hours)
#==========================================
#!/bin/bash
#SBATCH --job-name=blocks_parents
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=blocks_%A_%a.log
#SBATCH --array=0-8

set -euo pipefail

# ========= MODULES =========
module load StdEnv/2023
module load samtools/1.20
module load bedtools/2.31.0
module load picard
module load bamtools
module load java/17.0.6

# ========= USER VARIABLES =========
BASE_DIR="/home/yhezema/scratch/MAGIC"
ASSEMBLY_DIR="$BASE_DIR/Assembly/assemblies"
LOG_DIR="$BASE_DIR/Assembly/logs/rb_denovo"
REF_GENOME="$BASE_DIR/reference/Gmax_880_v6.0.fa"
REF_FAI="$BASE_DIR/reference/Gmax_880_v6.0.fai"

# Array of parent IDs
parents=(1355 Prosper RG11 RG22 RG23 RG64 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

WORK="$ASSEMBLY_DIR/${parent}"
LOG="$LOG_DIR/${parent}_map.log"

mkdir -p "$LOG_DIR"

echo "[Step 2+] Starting filtering/coverage for $parent at $(date)" >> "$LOG"

# ========= STEP 2.4: Split mapped/unmapped =========
echo "[Step 2.4] $parent Splitting mapped/unmapped reads started at [$(date)]" >> "$LOG"
samtools view -b -F 4 "$WORK/${parent}_mapped_all.bam" > "$WORK/${parent}_mapped_only.bam" 2>> "$LOG"
samtools index "$WORK/${parent}_mapped_only.bam" 2>> "$LOG"

samtools view -b -f 4 "$WORK/${parent}_mapped_all.bam" > "$WORK/${parent}_unmapped.bam" 2>> "$LOG"
samtools index "$WORK/${parent}_unmapped.bam" 2>> "$LOG"

# Split unmapped into paired only
samtools view -b -f 9 "$WORK/${parent}_unmapped.bam" > "$WORK/${parent}_unmapped_pair.bam" 2>> "$LOG"

# ========= STEP 2.5: Export unmapped reads in fastq =========
echo "[Step 2.5] Exporting unmapped reads to FASTQ at $(date)" >> "$LOG"

# Export properly paired unmapped reads using Picard
java -Xmx8G -jar $EBROOTPICARD/picard.jar SamToFastq \
  I="$WORK/${parent}_unmapped_pair.bam" \
  FASTQ="$WORK/${parent}_unmapped_R1.fastq" \
  SECOND_END_FASTQ="$WORK/${parent}_unmapped_R2.fastq" \
  VALIDATION_STRINGENCY=SILENT 2>> "$LOG"
rm "$WORK/${parent}_unmapped_pair.bam"

# Export unpaired unmapped reads using samtools
samtools fastq -1 "$WORK/${parent}_unmapped_unpair_R1.fastq" \
               -2 "$WORK/${parent}_unmapped_unpair_R2.fastq" \
               -0 /dev/null \
               -s /dev/null \
               -N \
               "$WORK/${parent}_unmapped.bam" 2>> "$LOG"

# ========= STEP 2.6: Filter mapped reads by mapping quality =========
echo "[Step 2.6] $parent Filtering mapped reads by MAPQ >= 10 started at $(date)" >> "$LOG"
samtools view -b -q 10 "$WORK/${parent}_mapped_only.bam" > "$WORK/${parent}_mapped_filtered.bam" 2>> "$LOG"
samtools index "$WORK/${parent}_mapped_filtered.bam" 2>> "$LOG"

# ========= STEP 2.7: Coverage for properly paired reads =========
echo "[Step 2.7] Generating coverage bedgraph (properly paired only) at $(date)" >> "$LOG"
samtools view -bf 0x2 "$WORK/${parent}_mapped_filtered.bam" \
  | samtools sort -n -@4 \
  | bedtools bamtobed -i - -bedpe \
  | awk '$1==$4' | cut -f1,2,6 | sort -k1,1 \
  | bedtools genomecov -i - -bga -g "${REF_GENOME}".fai > "$WORK/${parent}_coverage_paired.bedgraph" 2>> "$LOG"

# ========= STEP 2.8: Detect blocks =========
echo "[Step 2.8] $parent Detecting blocks with coverage >= 10 and merge distance 300 started at $(date)" >> "$LOG"
awk '$4 >= 10' "$WORK/${parent}_coverage_paired.bedgraph" > "$WORK/${parent}_coverage_filtered.bedgraph" 2>> "$LOG"
bedtools merge -i "$WORK/${parent}_coverage_filtered.bedgraph" -d 300 -c 4 -o mean > "$WORK/${parent}_blocks.txt" 2>> "$LOG"

# ========= STEP 2.9: Create superblocks =========
echo "[Step 2.10] $parent Creating superblocks (>= 12kb) started at $(date)" >> "$LOG"

# First get all blocks that are already >= 12kb
awk '{if ($3-$2 >= 12000) print $0}' "$WORK/${parent}_blocks.txt" > "$WORK/${parent}_superblocks.txt" 2>> "$LOG"

# Then process smaller blocks by merging them
bedtools merge -i "$WORK/${parent}_blocks.txt" -d 300 -c 4 -o mean \
  | awk '{if ($3-$2 >= 12000) print $0}' >> "$WORK/${parent}_superblocks.txt" 2>> "$LOG"

# Sort the final superblocks by chromosome and position
sort -k1,1 -k2,2n "$WORK/${parent}_superblocks.txt" > "$WORK/${parent}_superblocks_sorted.txt" 2>> "$LOG"

mv "$WORK/${parent}_superblocks_sorted.txt" "$WORK/${parent}_superblocks.txt"

# Remove any potential duplicates that might have been created
awk '!seen[$1,$2,$3]++' "$WORK/${parent}_superblocks.txt" > "$WORK/${parent}_superblocks_dedup.txt" 2>> "$LOG"

mv "$WORK/${parent}_superblocks_dedup.txt" "$WORK/${parent}_superblocks.txt"

echo "[Step 2.9] $parent Created $(wc -l < "$WORK/${parent}_superblocks.txt") superblocks at $(date)" >> "$LOG"









