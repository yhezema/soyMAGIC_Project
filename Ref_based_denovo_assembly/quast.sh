#!/bin/bash
#SBATCH --job-name=QUAST_10k
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00


### Load required modules (customize for your cluster)
module --force purge
module load python
#module load quast
#module spider quast
module load StdEnv/2020 gcc/9.3.0 python/3.8
module load quast/5.2.0

### Variables
PARENT="RG23"
THREADS=$SLURM_CPUS_PER_TASK
WORK="/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction
SCAFFOLDS="${WORK}/gapclosing/scaffolds_10k.fa
REF="/home/yhezema/scratch/MAGIC/reference/Gmax_880_v6.0.fa"
OUTPUT="${WORK}/final_quality/quast_10k"
LOG="${OUTPUT}/quast_${SLURM_JOB_ID}.log"

### Initialize
mkdir -p "$OUTPUT"
echo "=== QUAST JOB STARTED $(date) ===" > "$LOG"
echo "Job ID: $SLURM_JOB_ID" >> "$LOG"
echo "Assembly: $SCAFFOLDS" >> "$LOG"
echo "Reference: $REF" >> "$LOG"

### Run QUAST with optimized parameters
echo -e "\nRunning QUAST (minimap2 alignments)..." >> "$LOG"
quast.py \
  -o "$OUTPUT" \
  -r "$REF" \
  --eukaryote \
  --large \
  --min-alignment 500 \
  --no-icarus \
  --gene-finding \
  --fragmented \
  --labels "$PARENT" \
  "$SCAFFOLDS" \
  2>&1 | tee -a "$LOG"

### Check exit status
if [ $? -eq 0 ]; then
  echo -e "\nQUAST completed successfully at $(date)" >> "$LOG"
  echo "Results saved to: $OUTPUT" >> "$LOG"

  # Extract key metrics
  echo -e "\n=== KEY METRICS ===" >> "$LOG"
  grep -E "N50|L50|Total length|Largest contig|# N's per 100 kbp" "${OUTPUT}/report.txt" >> "$LOG"
else
  echo -e "\nERROR: QUAST failed (see ${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err)" >> "$LOG"

  # Fallback to reference-free mode if alignment fails
  echo -e "\nAttempting reference-free analysis..." >> "$LOG"
  quast.py \
    -o "${OUTPUT}_no_ref" \
    -t "$THREADS" \
    --eukaryote \
    --large \
    "$SCAFFOLDS" \
    2>&1 | tee -a "$LOG"
fi

### Finalize
echo -e "\n=== JOB COMPLETED $(date) ===" >> "$LOG"
#======================================================================
#======================================================================
