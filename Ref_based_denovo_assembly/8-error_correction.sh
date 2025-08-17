#======================================
# Step 8: error correction - Mapping reads to supercontigs
#======================================
# 1- Map reads to the indexed merged_supercontigs 
# 2- Error correction using racon and merged supercontigs as a reference and bam reads out from step #1
#==================================================
# 1. Mapping reads to supercontigs
#======================================

#!/bin/bash
#SBATCH --job-name=mapping_correct_map
#SBATCH --time=16:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=mapping_corr_%j.log
#SBATCH --error=mapping_corr_%j.err


# ===== Modules load =====
module load StdEnv/2023
module load bowtie2/2.5.4
module load samtools
module load picard

# ===== SETUP =====
parent="RG23"
WORK="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction"
TRIM_DIR="/home/yhezema/scratch/MAGIC/Assembly/trimmed/${parent}"
MERGED="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction/merged_supercontigs_500.fa"
REF_IDX_DIR="${WORK}/bowtie2_index_step8"
REF_IDX="${REF_IDX_DIR}/ref_index"
LOG_DIR="${WORK}/logs"
mkdir -p "$LOG_DIR"

ANALYSIS_LOG="$LOG_DIR/analysis.log"
ERROR_LOG="$LOG_DIR/error.log"

MAPPED_DIR=$WORK/mapping2
mkdir -p "$MAPPED_DIR"

sam_all="${MAPPED_DIR}/mapped_all.sam"
bam_all="${MAPPED_DIR}/mapped_all.bam"
bam_mapped="${MAPPED_DIR}/mapped.bam"
bam_filtered="${MAPPED_DIR}/mapped_filtered.bam"


R1_pe="${TRIM_DIR}/${parent}_R1_trimmed_pe.fastq.gz"
R2_pe="${TRIM_DIR}/${parent}_R2_trimmed_pe.fastq.gz"

cd "$WORK"

cp "/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/merged_supercontigs.fa" "$WORK"

# Filter merged superconyigs
seqkit seq -m 500 $WORK/merged_supercontigs.fa > "$MERGED"
seqkit stat $MERGED >> "$ANALYSIS_LOG"

echo "===== ERROR CORRECTION STARTED for ${parent} =====" | tee -a "$ANALYSIS_LOG"
date | tee -a "$ANALYSIS_LOG"
echo "Job ID: $SLURM_JOB_ID" | tee -a "$ANALYSIS_LOG"
echo "Host: $(hostname)" | tee -a "$ANALYSIS_LOG"

# ===== Build Bowtie2 index =====

echo "[$(date)] Building Bowtie2 index..."
bowtie2-build "$MERGED" "$REF_IDX" --threads 8 2>> "$ERROR_LOG"


# ===== 1. Map reads with Bowtie2 (SAM output) =====
echo "[$(date)] Mapping reads with Bowtie2 --sensitive..." | tee -a "$ANALYSIS_LOG"

bowtie2 --very-sensitive -p 16 \
    -q --phred33 \
    -x "$REF_IDX" \
    -1 "$R1_pe" \
    -2 "$R2_pe" \
    -S "$sam_all" \
    2>> "$ERROR_LOG"

if [[ ! -s "$sam_all" ]]; then
    echo "ERROR: SAM file not created or empty!" | tee -a "$ERROR_LOG"
    exit 1
fi

# ===== 2. Convert SAM to BAM and sort =====
echo "[$(date)] Converting SAM to sorted BAM..." | tee -a "$ANALYSIS_LOG"

samtools view -@ 8 -bS "$sam_all" \
    | samtools sort -@ 8 -T "${bam_all%.bam}_temp" -o "$bam_all" 2>> "$ERROR_LOG"

if [[ ! -s "$bam_all" ]]; then
    echo "ERROR: BAM file not created or empty!" | tee -a "$ERROR_LOG"
    exit 1
fi

samtools index "$bam_all"

# Remove sam file

rm "$sam_all"

# ===== 3. Mapping statistics =====
echo "All mapped reads statistics:" >> "$ANALYSIS_LOG"
samtools flagstat "$bam_all" >> "$ANALYSIS_LOG"

# ===== 4. Filter mapped reads only =====
samtools view -b -F 4 "$bam_all" > "$bam_mapped"
samtools index "$bam_mapped"
echo "Mapped reads statistics:" >> "$ANALYSIS_LOG"
samtools flagstat "$bam_mapped" >> "$ANALYSIS_LOG"

# ===== 5. Filter by quality >=10 =====
samtools view -b -q 10 "$bam_mapped" > "$bam_filtered"
samtools index "$bam_filtered"
echo "Quality filtered reads statistics:" >> "$ANALYSIS_LOG"
samtools flagstat "$bam_filtered" >> "$ANALYSIS_LOG"

echo "===== ERROR MAPPING COMPLETED for ${parent} =====" | tee -a "$ANALYSIS_LOG"
date | tee -a "$ANALYSIS_LOG"

#=========================================================================
#=========================================================================
# Extract fq reads
#-----------------
#!/bin/bash
#SBATCH --job-name=modern_error_correct
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=modern_corr_%j.log
#SBATCH --error=modern_corr_%j.err

set -Eeuo pipefail

# ===== Modules load =====
module load StdEnv/2023
module load samtools
module load python
module load bedtools
module load seqkit

module load racon 2>/dev/null 
module load minimap2 2>/dev/null 
module load bwa 2>/dev/null 

# ===== SETUP =====
parent="RG23"
WORK="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction"
CORRECTION="$WORK/correction"; mkdir -p "$CORRECTION"

# Use the correct original merged assembly as reference
MERGED="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction/merged_supercontigs_500.fa"

LOG_DIR="${WORK}/logs"; mkdir -p "$LOG_DIR"
EXEC_LOG="$LOG_DIR/execution.log"
ANALYSIS_LOG="$LOG_DIR/analysis.log"
ERROR_LOG="$LOG_DIR/error.log"

MAPPED_DIR="$WORK/mapping2"; mkdir -p "$MAPPED_DIR"
bam_filtered="${MAPPED_DIR}/mapped_filtered.bam"  # from previous job

threads="${SLURM_CPUS_PER_TASK:-8}"

log_msg "Starting modern error correction pipeline for $parent"
log_msg "Using reference: $MERGED"
log_msg "Threads: $threads"

# ===== INPUT VALIDATION =====
log_msg "Validating input files..."

check_file "$MERGED" "Reference assembly validation"
check_file "$bam_filtered" "Filtered BAM validation"

# ===== 1. Extract reads from BAM =====
log_msg "Extracting high-quality reads from BAM..."

fq1="$MAPPED_DIR/mapped_filtered_R1.fastq"
fq2="$MAPPED_DIR/mapped_filtered_R2.fastq"

# Extract paired reads with high quality
samtools fastq -@ "$threads" -1 "$fq1" -2 "$fq2" -0 /dev/null -s /dev/null -n -F 3844 "$bam_filtered" 2>> "$ERROR_LOG"

check_file "$fq1" "FASTQ R1 extraction"
check_file "$fq2" "FASTQ R2 extraction"

# Log read counts
r1_reads=$(( $(wc -l < "$fq1") / 4 ))
{
    echo "===== READ EXTRACTION STATS ====="
    echo "Timestamp: $(date)"
    echo "Extracted paired reads: $r1_reads"
} >> "$ANALYSIS_LOG"
log_msg "Extracted $r1_reads paired reads"
#======================================================================
#======================================================================
# Error correction - racon
#==========================
#!/bin/bash
#SBATCH --job-name=racon_chunks
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=racon_chunk_%j.log
#SBATCH --error=racon_chunk_%j.err

module load racon
module load seqkit   # for splitting FASTQ
module load minimap2 # to generate PAF for each chunk

WORK="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction"
CORRECTION="$WORK/correction"
MERGED="$WORK/merged_supercontigs_500.fa"
READS="$WORK/mapping2/all_reads.fastq"
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$CORRECTION/chunks"

# 1. Split the reads into chunks (e.g., 5M reads per chunk)
seqkit split2 -p 10 "$READS" -O "$CORRECTION/chunks"

# 2. Process each chunk sequentially
CHUNK_OUTS=()
for chunk in "$CORRECTION/chunks"/*.fastq; do
    base=$(basename "$chunk" .fastq)
    paf="$CORRECTION/${base}.paf"

    # Align chunk reads to merged contigs
    minimap2 -x map-ont -t "$THREADS" "$MERGED" "$chunk" > "$paf"

    # Run Racon for this chunk
    out_chunk="$CORRECTION/corrected_${base}.fa"
    racon -t "$THREADS" "$chunk" "$paf" "$MERGED" > "$out_chunk"

    CHUNK_OUTS+=("$out_chunk")
done
#====================================================================
#====================================================================
# Polishing - racon iteration
#=============================
#!/bin/bash
#SBATCH --job-name=polish_chunks
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-10
#SBATCH --output=polish_%A_%a.out
#SBATCH --error=polish_%A_%a.err

module load racon
module load minimap2

# Paths
WORK="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23"
CORRECTION="$WORK/error_correction"
READS_CHUNKS="$CORRECTION/chunks"

# Get chunk number with leading zeros (e.g., 001, 002)
CHUNK_ID=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

# File names
CONTIG="${CORRECTION}/corrected_all_reads.part_${CHUNK_ID}.fa"
READS="${READS_CHUNKS}/reads_chunk_$(printf "%02d" $SLURM_ARRAY_TASK_ID)"

# Polishing rounds
for ROUND in 1 2 3; do
    echo "=== Round $ROUND for chunk $CHUNK_ID ==="

    # Output file names
    if [ $ROUND -eq 1 ]; then
        INPUT_CONTIG="$CONTIG"
    else
        INPUT_CONTIG="${CORRECTION}/corrected_round$((ROUND-1))_${CHUNK_ID}.fa"
    fi
    OUTPUT_CONTIG="${CORRECTION}/corrected_round${ROUND}_${CHUNK_ID}.fa"

    # Map reads to contigs
    minimap2 -t $SLURM_CPUS_PER_TASK -x map-ont "$INPUT_CONTIG" "$READS" > "$CORRECTION/mappings_round${ROUND}_${CHUNK_ID}.paf"

    # Run Racon
    racon -t $SLURM_CPUS_PER_TASK "$READS" "$CORRECTION/mappings_round${ROUND}_${CHUNK_ID}.paf" "$INPUT_CONTIG" > "$OUTPUT_CONTIG"

    echo "Output saved to $OUTPUT_CONTIG"
done

seqkit fx2tab corrected_round3_merged.fa | wc -l
#=====================================================================================
#chunk-merging 
#===============
# keeps only the highest-quality corrected sequence for each region without duplication.

# Merge round 3 

CORRECTION=/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction/polish
cat $(printf "$CORRECTION/corrected_round3_%03d.fa " {1..10}) > $CORRECTION/corrected_round3_merged.fa

# Remove exact duplicates
seqkit rmdup corrected_round3_merged.fa > corrected_round3_dedup.fa
#[INFO] 2260122 duplicated records removed
seqkit stats corrected_round3_dedup.fa
 
# Keep only sequences >= 500bp 
seqkit seq -m 200 corrected_round3_dedup.fa > corrected_round3_200.fa
seqkit stats corrected_round3_200.fa

# Keep only sequences >= 500bp 
seqkit seq -m 500 corrected_round3_dedup.fa > corrected_round3_500.fa
seqkit stats corrected_round3_500.fa

seqkit seq -m 1000 corrected_round3_dedup.fa > corrected_round3_1k.fa
seqkit stats corrected_round3_1k.fa

GAPCLOSING="$WORK/error_correction/gapclosing"
mkdir -p $GAPCLOSING

cd $CORRECTION/corrected_round3_merged.fa $GAPCLOSING
#====================================================================
#====================================================================

