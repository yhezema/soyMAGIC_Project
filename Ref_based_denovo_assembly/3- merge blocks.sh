#==========================================================
# Step 3: Defining blocks and superblocks (merge blocks) (4-5 hours)
#==========================================================
#!/bin/bash
#SBATCH --job-name=blocks
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=blocks_%A_%a.log
#SBATCH --array=0-7

# Array of parent IDs
parents=(1355 Prosper RG11 RG22 RG23 RG64 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# To ensures script stops on any error
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
  | bedtools genomecov -i - -bga -g "$REF_GENOME".fai > "$WORK/${parent}_coverage_paired.bedgraph" 2>> "$LOG"

# ========= STEP 2.8: Detect blocks =========
echo "[Step 2.8] $parent Detecting blocks with coverage >= 10 and merge distance 300 started at $(date)" >> "$LOG"
awk '$4 >= 10' "$WORK/${parent}_coverage.bedgraph" > "$WORK/${parent}_coverage_filtered.bedgraph" 2>> "$LOG"
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

#=============================================================================
