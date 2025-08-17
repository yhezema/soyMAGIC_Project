#====================================================================
------------------------------
# Step 9: Run RagTag for gap closing
-----------------------------
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=ragtag_final
#SBATCH --output=ragtag_final.log
#SBATCH --error=ragtag.err

# List of parents
parent=RG23

# Activate virtual environment
module load StdEnv/2023 minimap2/2.28

source ~/ragtag-venv/bin/activate

# Variables
REF=/home/yhezema/scratch/MAGIC/reference/Gmax_880_v6.0.fa
BASE=/home/yhezema/scratch/MAGIC/Assembly/assemblies/RG23/error_correction
CORRECTED=$BASE/polish/corrected_round3_200.fa
OUTDIR=$BASE/gapclosing
#THREADS=${SLURM_CPUS_PER_TASK:-8}
THREADS=16

mkdir -p $OUTDIR

# Run RagTag
ragtag.py scaffold -u "$REF" "$CORRECTED" -o "$OUTDIR" -t "$THREADS"

seqkit seq -m 1000  $OUTDIR/ragtag.scaffold.fasta > scaffolds_1k.fa
seqkit stats $OUTDIR/scaffolds_1k.fa

seqkit seq -m 10000  $OUTDIR/ragtag.scaffold.fasta > scaffolds_10k.fa
seqkit stats $OUTDIR/scaffolds_10k.fa

#====================================================================
#====================================================================
