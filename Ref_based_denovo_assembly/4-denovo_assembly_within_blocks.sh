#====================================
# Step 4: Assembly within superblocks
#====================================
#!/bin/bash
#SBATCH --job-name=assembly_mapped
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=Assembly_mapped_%A_%a.log
#SBATCH --error=Assembly_mapped_%j.err
#SBATCH --array=0-7

# Array of parent IDs
parents=(1355 Prosper RG11 RG22 RG23 RG64 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}


echo "===================================================================="
echo "SUPERBLOCK ASSEMBLY_RG46"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "===================================================================="

# ===========================================
# Assemble within  Superblocks  (12-16 h)
# ===========================================
# To ensures script stops on any error
set -euo pipefail

# ----------------------------
# Initialize Environment
# ----------------------------
BASE_DIR="/home/yhezema/scratch/MAGIC"
WORK="$BASE_DIR/Assembly/assemblies/${parent}"
ERROR_LOG="$WORK/${parent}_assembly.error.log"
LOG="$WORK/${parent}_assembly.log"

# Create output directory
mkdir -p "$WORK"/{abyss_results/{contigs,scaffolds,stats},unassembled_reads}
touch "$ERROR_LOG"

# ----------------------------
# Load Required Modules
# ----------------------------
module --force purge
module load StdEnv/2023 gcc/12.3 abyss/2.3.7 samtools bedtools bamtools seqtk

echo -e "\n===== STARTING FULL ASSEMBLY =====" | tee -a "$LOG"
echo "Start time: $(date)" | tee -a "$LOG"

# -------------------------------------
# Filter and split superblocks >100000
# -------------------------------------
{
  echo -e "\n===== [STEP 3.1] VALIDATING INPUTS AND FILTERING SUPERBLOCKS ====="
  echo "Start time: $(date)"

   awk '{if ($3-$2 >= 100000) {split_size=100000; overlap=300; for(start=$2; start<$3; start+=split_size-overlap) {end=(start+split_size>$3)?$3:start+split_size; print $1":"start"-"end}} else print $1":"$2"-"$3}' "$WORK/${parent}_superblocks.txt" > "$WORK/filtered_blocks.txt"


  echo "Filtered blocks saved to: $WORK/filtered_blocks.txt"
  echo "Step completed at: $(date)"
} | tee -a "$LOG" 

# ----------------------------
# Final process_block Function
# ----------------------------
process_block() {
  block=$1
  block_id=$(echo "$block" | tr ':' '_' | tr '-' '_')
  temp_dir=$(mktemp -d -p "$WORK")

 # Skip block if contigs, scaffolds, and stats already exist
  contig_file="$WORK/abyss_results/contigs/block_${block_id}-contigs.fa"
  scaffold_file="$WORK/abyss_results/scaffolds/block_${block_id}-scaffolds.fa"
  stats_file="$WORK/abyss_results/stats/block_${block_id}-stats.csv"

  if [[ -s "$contig_file" && -s "$scaffold_file" && -s "$stats_file" ]]; then
    echo "Block $block_id already assembled (contigs, scaffolds, stats present). Skipping..." | tee -a "$LOG"
    return 0
  fi

  echo "[$(date)] Processing block: $block_id" | tee -a "$LOG"

  # 1. Extract reads from sorted BAM
  echo "[DEBUG] Extracting reads for $block_id" | tee -a "$LOG"

  if ! samtools view -b "$WORK/${parent}_mapped_filtered.bam" "$block" > "$temp_dir/block.bam"; then
    echo "ERROR: Failed to extract reads by name" | tee -a "$ERROR_LOG"
    return 1
  fi

  # 2. Read statistics
  echo "[DEBUG] Read statistics for $block_id:" | tee -a "$LOG"
 
  total_reads=$(samtools view -c "$temp_dir/block.bam")
  proper_pairs=$(samtools view -c -f 2 "$temp_dir/block.bam")
  singleton=$(( total_reads - proper_pairs * 2 ))
  [[ $singleton -lt 0 ]] && singleton=0
  block_size=$(echo "$block" | awk -F"[:-]" '{print $3-$2}')
  echo "Total reads: $total_reads" | tee -a "$LOG"
  echo "Proper pairs: $proper_pairs" | tee -a "$LOG"
  echo "Singleton reads: $singleton" | tee -a "$LOG"
  echo "Block size: $block_size bp" | tee -a "$LOG"
  
  
  echo "Read density: $(awk -v r=$total_reads -v s=$block_size 'BEGIN{printf "%.1f", r*150/s}')x" | tee -a "$LOG"

  if [[ $total_reads -lt 1000 ]]; then
    echo "WARNING: Block $block_id has $total_reads reads (<1000)" | tee -a "$ERROR_LOG"
  fi

  # 3. Convert to FASTQ
  echo "[DEBUG] Converting to FASTQ" | tee -a "$LOG"
  bedtools bamtofastq -i "$temp_dir/block.bam" \
    -fq "$temp_dir/R1.fq" \
    -fq2 "$temp_dir/R2.fq" 2> /dev/null

  # 4. Recover singleton reads
 echo "[DEBUG] Recovering singletons" | tee -a "$LOG"
 {
  samtools view -b -f 72 -F 8 "$temp_dir/block.bam" \
    | bamtools convert -format fastq > "$temp_dir/singleton_R1.fq" 2>> "$ERROR_LOG" || touch "$temp_dir/singleton_R1.fq"

  samtools view -b -f 136 -F 8 "$temp_dir/block.bam" \
    | bamtools convert -format fastq > "$temp_dir/singleton_R2.fq" 2>> "$ERROR_LOG" || touch "$temp_dir/singleton_R2.fq"
 } || true

 # 5. Recover orphaned mates skipped by bedtools
 echo "[DEBUG] Recovering orphaned non-adjacent mates" | tee -a "$LOG"
 {
  samtools view -b -f 64 -F 2 "$temp_dir/block.bam" \
    | bamtools convert -format fastq > "$temp_dir/orphan_R1.fq" 2>> "$ERROR_LOG" || touch "$temp_dir/orphan_R1.fq"

  samtools view -b -f 128 -F 2 "$temp_dir/block.bam" \
    | bamtools convert -format fastq > "$temp_dir/orphan_R2.fq" 2>> "$ERROR_LOG" || touch "$temp_dir/orphan_R2.fq"
 } || true

 # 6. Combine all read types into R1 and R2
 echo "[DEBUG] Combining all reads" | tee -a "$LOG"
 cat "$temp_dir/R1.fq" "$temp_dir/singleton_R1.fq" "$temp_dir/orphan_R1.fq" > "$temp_dir/all_R1.fq"
 cat "$temp_dir/R2.fq" "$temp_dir/singleton_R2.fq" "$temp_dir/orphan_R2.fq" > "$temp_dir/all_R2.fq"

  # 7. Run ABYSS

  echo "[DEBUG] Running ABYSS for $block_id" | tee -a "$LOG"
  start_time=$(date +%s)

  # Sanity check: Ensure FASTQ files are not empty
  if [[ ! -s "$temp_dir/all_R1.fq" || ! -s "$temp_dir/all_R2.fq" ]]; then
    echo "WARNING: Skipping $block_id due to empty FASTQ files." | tee -a "$ERROR_LOG"
    return 1
  fi

  # Run ABYSS
  if abyss-pe k=63 B=4G np=4 j=4 s=100 \
    name="$temp_dir/block_${block_id}" \
    in="$temp_dir/all_R1.fq $temp_dir/all_R2.fq" 2>> "$ERROR_LOG"; then
  
    echo "[SUCCESS] Assembled $block_id" | tee -a "$LOG"

  # Copy contigs
   if [[ -f "$temp_dir/block_${block_id}-contigs.fa" ]]; then
     cp "$temp_dir/block_${block_id}-contigs.fa" "$WORK/abyss_results/contigs/" 2>> "$ERROR_LOG"
     contig_count=$(grep -c '^>' "$temp_dir/block_${block_id}-contigs.fa")
   else
     echo "WARNING: Contigs file missing for $block_id" | tee -a "$ERROR_LOG"
     contig_count=0
   fi

   # Copy scaffolds
   if [[ -f "$temp_dir/block_${block_id}-scaffolds.fa" ]]; then
     cp "$temp_dir/block_${block_id}-scaffolds.fa" "$WORK/abyss_results/scaffolds/" 2>> "$ERROR_LOG"
     scaffold_count=$(grep -c '^>' "$temp_dir/block_${block_id}-scaffolds.fa")
     gap_bases=$(grep -o 'N\+' "$temp_dir/block_${block_id}-scaffolds.fa" | awk '{sum += length($0)} END {print sum + 0}')

     echo "Gap bases in scaffolds: ${gap_bases}" | tee -a "$LOG"
   else
     echo "No scaffolds generated for $block_id" | tee -a "$LOG"
     scaffold_count=0
   fi

   # Copy stats
   if [[ -f "$temp_dir/block_${block_id}-stats.csv" ]]; then
     cp "$temp_dir/block_${block_id}-stats.csv" "$WORK/abyss_results/stats/" 2>> "$ERROR_LOG"
   fi

   # Log summary
   echo "Assembly stats for $block_id: $contig_count contigs, $scaffold_count scaffolds" | tee -a "$LOG"

 else
   echo "ERROR: ABYSS assembly failed for $block_id" | tee -a "$ERROR_LOG"
   return 1
 fi

 end_time=$(date +%s)
 runtime=$((end_time - start_time))
 echo "Runtime for $block_id: ${runtime}s" | tee -a "$LOG"

 # Cleanup
 rm -rf "$temp_dir"

}
export WORK
export parent
export LOG
export ERROR_LOG
export -f process_block

# ----------------------------
# Run parallel assembly
# ----------------------------
# Convert input block list to Unix format (prevents read issues)
dos2unix "$WORK/filtered_blocks.txt" 2>/dev/null

# Count total numbers of blocks
echo "Processing $(wc -l < $WORK/filtered_blocks.txt) superblocks"

# Read each block and process
parallel --jobs 4 --joblog "$WORK/parallel_block_assembly.log" process_block :::: "$WORK/filtered_blocks.txt"

#-----------------------------
# Verify Results
# ----------------------------
echo -e "\n===== ASEMBLY RESULTS =====" | tee -a "$LOG"
echo "Contigs generated:" | tee -a "$LOG"
ls -lh "$WORK/abyss_results/contigs/" | tee -a "$LOG"
echo -e "\nScaffolds generated:" | tee -a "$LOG"
ls -lh "$WORK/abyss_results/scaffolds/" | tee -a "$LOG"
echo -e "\nStats files:" | tee -a "$LOG"
ls -lh "$WORK/abyss_results/stats/" | tee -a "$LOG"

echo -e "\n===== ASSEMBLY COMPLETED =====" | tee -a "$LOG"
echo "End time: $(date)" | tee -a "$LOG"
echo "Check $LOG and $ERROR_LOG for details" | tee -a "$LOG"
#===========================================================================================================
#===========================================================================================================
# =============================================================
# Unmapped Reads Assembly 
# =============================================================
#!/bin/bash
#SBATCH --job-name=assem_unmapped
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=assem_unmap_%A_%a.log
#SBATCH --error=assem_unmap_%j.err
#SBATCH --array=0-7

# To ensures script stops on any error
set -euo pipefail

# List of parents
parents=(1355 Prosper RG11 RG22 RG23 RG46 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# ----------------------------
# Variables
# ----------------------------
BASE_DIR="/home/yhezema/scratch/MAGIC"
WORK="$BASE_DIR/Assembly/assemblies/${parent}"
LOG_DIR="$WORK/logs"
ERROR_LOG="$LOG_DIR/${parent}_assembly.error.log"
LOG="$LOG_DIR/${parent}_assembly.log"
UNMAPPED_DIR="$WORK/unassembled_reads"
RESULTS_DIR="$WORK/abyss_results/unmap_contigs"
FINAL_CONTIGS="$RESULTS_DIR/${parent}_final_contigs.fa"

# Create directories
mkdir -p "$UNMAPPED_DIR" "$RESULTS_DIR" "$LOG_DIR"
touch "$ERROR_LOG" "$LOG"

# ----------------------------
# Load Modules
# ----------------------------
module --force purge
module load StdEnv/2023 gcc/12.3 abyss seqtk seqkit

# ----------------------------
# Start the assembly
# ----------------------------
{
echo "===== UNMAPPED READS ASSEMBLY FOR ${parent} ====="
echo "Start: $(date)"
echo "Working dir: $WORK"
[[ -n "${SLURM_JOB_ID:-}" ]] && echo "Job ID: $SLURM_JOB_ID"

# Combine all reads
echo -e "\n[1] Combining reads..."
cat "$WORK/${parent}_unmapped_R1.fastq" \
    "$WORK/${parent}_unmapped_R2.fastq" \
    "$WORK/${parent}_unmapped_unpair_R1.fastq" \
    "$WORK/${parent}_unmapped_unpair_R2.fastq" \
    > "$UNMAPPED_DIR/${parent}_combined_reads.fastq"

# Run ABySS assembly
echo -e "\n[2] Running ABySS..."

abyss-pe k=63 c=1 name="$UNMAPPED_DIR/assembly" \
    in="$WORK/${parent}_unmapped_R1.fastq $WORK/${parent}_unmapped_R2.fastq" \
    se="$WORK/${parent}_unmapped_unpair_R1.fastq $WORK/${parent}_unmapped_unpair_R2.fastq" \
    j=4 B=16G np=4 s=100

# Process results
echo -e "\n[3] Processing contigs..."
if [[ -f "$UNMAPPED_DIR/assembly-contigs.fa" ]]; then
    cp "$UNMAPPED_DIR/assembly-contigs.fa" "$RESULTS_DIR/${parent}_contigs.fa"
    seqkit seq -m 500 "$RESULTS_DIR/${parent}_contigs_500.fa" > "$FINAL_CONTIGS"
else
    echo "ERROR: No contigs generated" | tee -a "$ERROR_LOG"
    exit 1
fi

# Final stats
echo -e "\n[4] Assembly statistics:"
seqkit stat "$FINAL_CONTIGS"
echo -e "\nCompleted: $(date)"
} | tee -a "$LOG"

# Final check
if [[ $? -eq 0 ]] && [[ -s "$FINAL_CONTIGS" ]]; then
    echo "SUCCESS: Assembly completed" | tee -a "$LOG"
    exit 0
else
    echo "ERROR: Assembly failed" | tee -a "$ERROR_LOG"
    exit 1
fi
#=============================================================================

