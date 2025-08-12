#=================================================
# Step 6: Mapping to Supercontigs  (12-16 hours)
#=================================================
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=map_supcont
#SBATCH --output=map_supcont_%A_%a.log
#SBATCH --error=map_supcont.err
#SBATCH --array=0-7

# Array of parent IDs
parents=(1355 Prosper RG11 RG22 RG23 RG64 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# ========= MODULES =========
module load StdEnv/2023
module load samtools/1.20
module load bedtools/2.31.0
module load picard
module load bamtools
module load java/17.0.6
module load bowtie2/2.5.4
module load abyss
module load seqkit

# Initialize variables
BASE="/home/yhezema/scratch/MAGIC/Assembly/assemblies/${parent}"
TRIM_DIR="/home/yhezema/scratch/MAGIC/Assembly/trimmed/${parent}"

# Initialize logging
LOG_DIR="${BASE}/logs"
mkdir -p "$LOG_DIR"
EXEC_LOG="$LOG_DIR/step5_execution.log"
ERROR_LOG="$LOG_DIR/step5_error.log"
ANALYSIS_LOG="$LOG_DIR/step5_analysis.log"

# Move the scaffolds file to the base dir
cp "${BASE}/ragtag_output_5/ragtag.scaffold.fasta" "${BASE}/ragtag.scaffold_500.fasta"

{
echo "===== SUPERCONTIG RECOVERY PIPELINE STARTED for ${parent}  ====="
date
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"

# Define paths and filenames
SC_DIR="${BASE}/step5"
SUPERCONTIGS="${BASE}/ragtag.scaffold_500.fasta"
UNIQUE_SC="${BASE}/ragtag.scaffold_500.renamed.fasta"
REF_IDX_DIR="${BASE}/bowtie2_index"
REF_IDX="${REF_IDX_DIR}/ref_index"

UNMAPPED_DIR="${SC_DIR}/UnmappedReads"
ASSEMBLY_DIR="${SC_DIR}/Unassembled"
mkdir -p "$UNMAPPED_DIR" "$ASSEMBLY_DIR" "$(dirname "$REF_IDX")"

# Input files
R1_pe="${TRIM_DIR}/${parent}_R1_trimmed_pe.fastq.gz"
R2_pe="${TRIM_DIR}/${parent}_R2_trimmed_pe.fastq.gz"
R1_se="${TRIM_DIR}/${parent}_R1_trimmed_se.fastq.gz"
R2_se="${TRIM_DIR}/${parent}_R2_trimmed_se.fastq.gz"

# Validate input files
echo "===== INPUT VALIDATION =====" >> "$EXEC_LOG"
for file in "$SUPERCONTIGS" "$R1_pe" "$R2_pe" "$R1_se" "$R2_se"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing input file $file" >> "$ERROR_LOG"
        exit 1
    else
        echo "Found: $file" >> "$EXEC_LOG"
    fi
done

# Step 1: Make sequence names unique
echo "===== RENAMING CONTIGS =====" >> "$EXEC_LOG"
seqkit rename "$SUPERCONTIGS" > "$UNIQUE_SC" 2>> "$ERROR_LOG"
seqkit stat "$UNIQUE_SC" >> "$ANALYSIS_LOG"

# Step 2: Build Bowtie2 index
echo "===== BUILDING BOWTIE2 INDEX =====" >> "$EXEC_LOG"
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" "$UNIQUE_SC" "$REF_IDX" 2>> "$ERROR_LOG"

# Verify index creation
if [[ ! -f "${REF_IDX}.1.bt2" ]]; then
    echo "ERROR: Bowtie2 index files not created" >> "$ERROR_LOG"
    exit 1
fi

# Step 3: Map reads to supercontigs
echo "===== MAPPING READS =====" >> "$EXEC_LOG"
bowtie2 --very-sensitive \
    -p "$SLURM_CPUS_PER_TASK" \
    -x "$REF_IDX" \
    -1 "$R1_pe" -2 "$R2_pe" \
    -U "$R1_se,$R2_se" \
    -S "$UNMAPPED_DIR/mapped.sam" 2> "$UNMAPPED_DIR/mapping_stats.txt"

# Log mapping statistics
cat "$UNMAPPED_DIR/mapping_stats.txt" >> "$ANALYSIS_LOG"

# Step 4: Process unmapped reads
echo "===== PROCESSING UNMAPPED READS =====" >> "$EXEC_LOG"
{
    samtools view -@ "$SLURM_CPUS_PER_TASK" -bS "$UNMAPPED_DIR/mapped.sam" | \
    samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$UNMAPPED_DIR/mapped.bam"
    samtools index "$UNMAPPED_DIR/mapped.bam"

    # Extract unmapped reads
    samtools view -@ "$SLURM_CPUS_PER_TASK" -b -f 12 -F 256 "$UNMAPPED_DIR/mapped.bam" > "$UNMAPPED_DIR/unmapped_paired.bam"
    samtools view -@ "$SLURM_CPUS_PER_TASK" -b -f 4 "$UNMAPPED_DIR/mapped.bam" > "$UNMAPPED_DIR/unmapped_unpaired.bam"

    # Convert to FASTQ
    java -jar $EBROOTPICARD/picard.jar SamToFastq \
        I="$UNMAPPED_DIR/unmapped_paired.bam" \
        FASTQ="$UNMAPPED_DIR/unmapped_paired_R1.fastq" \
        SECOND_END_FASTQ="$UNMAPPED_DIR/unmapped_paired_R2.fastq" \
        VALIDATION_STRINGENCY=SILENT 2>> "$ERROR_LOG"

    samtools fastq -@ "$SLURM_CPUS_PER_TASK" "$UNMAPPED_DIR/unmapped_unpaired.bam" > "$UNMAPPED_DIR/unmapped_unpaired.fastq" 2>
} >> "$EXEC_LOG" 2>> "$ERROR_LOG"

# Uunmapped read counts
echo "===== UNMAPPED READ COUNTS =====" >> "$ANALYSIS_LOG"
echo "Paired unmapped: $(grep -c '^@' "$UNMAPPED_DIR/unmapped_paired_R1.fastq") pairs" >> "$ANALYSIS_LOG"
echo "Unpaired unmapped: $(grep -c '^@' "$UNMAPPED_DIR/unmapped_unpaired.fastq") reads" >> "$ANALYSIS_LOG"
} > >(tee -a "$EXEC_LOG") 2> >(tee -a "$ERROR_LOG" >&2)
#=====================================================================================================
