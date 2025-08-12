#============================================
# Step 5: Merge contigs and remove redundancy
#============================================
#!/bin/bash
#SBATCH --job-name=merge_contigs
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=merge_contigs_%A_%a.log
#SBATCH --error=merge_contigs_%j.err
#SBATCH --array=0-7

# List of parents
parents=(1355 Prosper RG11 RG22 RG23 RG46 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}

# To ensures script stops on any error
set -euo pipefail

# ----------------------------
# Environment and Variables
# ----------------------------
NThreads=8
BASE="/home/yhezema/scratch/MAGIC/Assembly/assemblies/${parent}"
REF="/home/yhezema/scratch/MAGIC/reference/Gmax_880_v6.0.fa"
ABYSS_DIR="$BASE/abyss_results"
UNASSEMBLED="$BASE/unassembled_reads"
RESULTS_DIR="$ABYSS_DIR/contigs"
WORK="$BASE"
FINAL_CONTIGS="$BASE/${parent}_final_contigs.fa"

LOG="$BASE/logs"
ERROR_LOG="$BASE/logs"
mkdir -p "$LOG" "$ERROR_LOG"

ERROR_LOG="$LOG_DIR/step4_supercontigs.error.log"
LOG="$LOG_DIR/step4_supercontigs.log"

# ----------------------------
# Load Modules (Cedar-specific)
# ----------------------------
module --force purge
module load StdEnv/2023
module load python/3.10.2
module load java
module load seqkit

# ----------------------------
# Step 5.1: Merge all contigs with QC
# ----------------------------

echo -e "\n[STEP 4.1] Merging and QC of contigs..." | tee -a "$LOG"


# Move unmapped contigs 
cp "$ABYSS_DIR/unmap_contigs/${parent}_contigs.fa" "$RESULTS_DIR/${parent}_contigs.fa"

# Merge all contigs into one file
cat "$RESULTS_DIR"/*.fa > "$WORK/combined_contigs.fa" 


    # Step 1: Rename, clean, and deduplicate
 seqkit rename "$WORK/combined_contigs.fa" | \
 seqkit seq -g | \
 seqkit rmdup -s | \
awk '
    BEGIN { OFS="\n"; i=0 }
    /^>/ { print ">contig_" ++i; next }
    {
      seq = toupper($0)
      gsub(/[^ACGTN]/, "N", seq)
      print seq
    }
' | \
seqkit sort -l -r > "$WORK/renamed_cleaned_contigs.fa"

# Filter final contigs
seqkit seq -m 500 "$WORK/renamed_cleaned_contigs.fa" > "$FINAL_CONTIGS" 

# Verify output
echo "Final contig count: $(grep -c '^>' "$FINAL_CONTIGS")" 

---
# Verify
seqkit stat "$WORK/renamed_cleaned_contigs.fa" >> "$LOG"
seqkit stat "$WORK/${parent}_final_contigs.fa" >> "$LOG"

============================================================
============================================================
------------------------------
# Step 5.2: Run RagTag to get supercontigs
-----------------------------
#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=ragtag
#SBATCH --output=ragtag_%A_%a.log
#SBATCH --error=ragtag.err

#SBATCH --array=0-7

# List of parents
parents=(1355 Prosper RG11 RG22 RG23 RG46 X790 07-78)
parent=${parents[$SLURM_ARRAY_TASK_ID]}


# Activate virtual environment
module load StdEnv/2023 minimap2/2.28

source ~/ragtag-venv/bin/activate

# Variables
REF=/home/yhezema/scratch/MAGIC/reference/Gmax_880_v6.0.fa
CONTIGS=/home/yhezema/scratch/MAGIC/Assembly/assemblies/${parent}/${parent}_final_contigs.fa
OUTDIR=/home/yhezema/scratch/MAGIC/Assembly/assemblies/${parent}/ragtag_output_5
#THREADS=${SLURM_CPUS_PER_TASK:-8}
THREADS=16

# Run RagTag
ragtag.py scaffold "$REF" "$CONTIGS" -o "$OUTDIR" -t "$THREADS"
    
#=============================================================================

