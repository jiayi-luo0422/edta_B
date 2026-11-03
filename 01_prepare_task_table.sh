#!/bin/bash
#SBATCH --job-name=edta_tasks
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/edta_tasks_%j.out
#SBATCH --error=logs/edta_tasks_%j.err

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/workflow_common.sh"
activate_edta
mkdir -p "$PROJDIR" "$(dirname "$TASK_TABLE")" "$(dirname "$SAMPLE_TABLE")"
samples_tmp=$(mktemp "$SAMPLE_TABLE.XXXXXX")
tasks_tmp=$(mktemp "$TASK_TABLE.XXXXXX")
trap 'rm -f "$samples_tmp" "$tasks_tmp"' EXIT
for genome in "$GENOME_DIR"/*.chrom.fasta; do
    [ -s "$genome" ] || continue
    basename "$genome" .chrom.fasta
done | LC_ALL=C sort > "$samples_tmp"
[ -s "$samples_tmp" ] || { echo "No genomes in $GENOME_DIR" >&2; exit 1; }
while IFS= read -r sample; do
    genome="$GENOME_DIR/$sample.chrom.fasta"
    if [ ! -s "$genome.fai" ] || [ "$genome" -nt "$genome.fai" ]; then
        "$SAMTOOLS_BIN" faidx "$genome"
    fi
    awk -F '\t' -v sample="$sample" 'NF>=2 {print sample "\t" $1; n++} END {if(!n) exit 1}' "$genome.fai" >> "$tasks_tmp"
done < "$samples_tmp"
# Keep mtimes when the prepared inputs have not changed.
if ! cmp -s "$samples_tmp" "$SAMPLE_TABLE"; then mv "$samples_tmp" "$SAMPLE_TABLE"; fi
if ! cmp -s "$tasks_tmp" "$TASK_TABLE"; then mv "$tasks_tmp" "$TASK_TABLE"; fi
printf 'Samples: %s\nChromosome tasks: %s\n' "$(wc -l < "$SAMPLE_TABLE")" "$(wc -l < "$TASK_TABLE")"
