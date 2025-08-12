=============================================
# Step 7: De Novo Assembly of Unmapped Reads 
=============================================
#!/bin/bash
#SBATCH --job-name=assemb_unmapped
#SBATCH --output=assemb_unmapped.out
#SBATCH --error=assemb_unmapped.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-7

# Array of parent IDs
parents=(1355 Prosper RG11 RG22 RG23 RG64 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# To avoid "Runtime.totalMemory()" error
export _JAVA_OPTIONS="-Xmx16g"

# ======= MODULES =======
module load StdEnv/2023
module load abyss/2.3.7
module load seqkit
module load samtools/1.20

# ======= VARIABLES =======
BASE="/home/yhezema/scratch/MAGIC/Assembly/assemblies/${parent}"
SC_DIR="${BASE}/step5"
UNMAPPED_DIR="${SC_DIR}/UnmappedReads"
ASSEMBLY_DIR="${SC_DIR}/Unassembled"
MERGED_CONTIGS="${ASSEMBLY_DIR}/${parent}_unmapped_merged.fa"

mkdir -p "$ASSEMBLY_DIR"

LOG_DIR="${BASE}/logs"
mkdir -p "$LOG_DIR"
EXEC_LOG="$LOG_DIR/step5_execution.log"
ERROR_LOG="$LOG_DIR/step5_error.log"
ANALYSIS_LOG="$LOG_DIR/step5_analysis.log"

# ======= STEP 5.0: Filter very short reads to save memory =======
echo "[$(date)] Filtering unmapped reads shorter than 50 bp..."
seqkit seq -m 50 "${UNMAPPED_DIR}/unmapped_paired_R1.fastq" > "${UNMAPPED_DIR}/unmapped_paired_R1.filt.fastq"
seqkit seq -m 50 "${UNMAPPED_DIR}/unmapped_paired_R2.fastq" > "${UNMAPPED_DIR}/unmapped_paired_R2.filt.fastq"
seqkit seq -m 50 "${UNMAPPED_DIR}/unmapped_unpaired.fastq"  > "${UNMAPPED_DIR}/unmapped_unpaired.filt.fastq"

# ======= STEP 5.1: Run ABySS for ${parent} =======
{
cd "$ASSEMBLY_DIR"

echo "[$(date)] Running ABySS assembly at k=63..."

abyss-pe name=Unmapped_k63 k=63 j="$SLURM_CPUS_PER_TASK" B=16G \
    lib="pe" \
    pe="$UNMAPPED_DIR/unmapped_paired_R1.filt.fastq $UNMAPPED_DIR/unmapped_paired_R2.filt.fastq" \
    se="$UNMAPPED_DIR/unmapped_unpaired.filt.fastq" \
    2>> "$ERROR_LOG"

# Assembly statistics
echo "===== ASSEMBLY STATISTICS =====" >> "$ANALYSIS_LOG"
seqkit stat "Unmapped_k63-contigs.fa" >> "$ANALYSIS_LOG"

# check unictigs stats
seqkit stat "$ASSEMBLY_DIR/Unmapped_k63-unitigs.fa" >> "$ANALYSIS_LOG"


# Filter contigs
echo "===== FILTERING CONTIGS =====" >> "$EXEC_LOG"
#seqkit seq -m 200 "$ASSEMBLY_DIR/Unmapped_k63-unitigs.fa" > "$ASSEMBLY_DIR/Unmapped_contigs_200.fa" 2>> "$ERROR_LOG"

#seqkit stat "$ASSEMBLY_DIR/Unmapped_contigs_200.fa" >> "$ANALYSIS_LOG"

seqkit seq -m 500 "$ASSEMBLY_DIR/Unmapped_k63-unitigs.fa" > "$ASSEMBLY_DIR/Unmapped_contigs_500.fa"

seqkit stat "$ASSEMBLY_DIR/Unmapped_contigs_500.fa" >> "$ANALYSIS_LOG"

echo "===== PIPELINE COMPLETED at $(date) ====="

#=============================================================================
