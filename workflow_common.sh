#!/bin/bash
WORKFLOW_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || return 1
source "${EDTA_CONFIG:-$WORKFLOW_DIR/edta_config.sh}" || return 1
# A user-supplied config may override only the fields it needs.
PROJDIR=${PROJDIR:?Set PROJDIR in EDTA_CONFIG}
GENOME_DIR=${GENOME_DIR:-$PROJDIR/data/chromosome_only}
EDTA_OUT_DIR=${EDTA_OUT_DIR:-$PROJDIR/data/edta_new}
CDS_DIR=${CDS_DIR:-$PROJDIR/cds}
PANTE_DIR=${PANTE_DIR:-$PROJDIR/04_panTE}
PANTE_MASK_DIR=${PANTE_MASK_DIR:-$PROJDIR/05_panTE_annotation}
TASK_TABLE=${TASK_TABLE:-$EDTA_OUT_DIR/edta_tasks.tsv}
SAMPLE_TABLE=${SAMPLE_TABLE:-$EDTA_OUT_DIR/edta_samples.txt}
CONDA_SH=${CONDA_SH:-/public/home/a326/luojy/miniconda3/etc/profile.d/conda.sh}
SAMTOOLS_BIN=${SAMTOOLS_BIN:-samtools}
EDTA_ENV=${EDTA_ENV:-EDTA2_fix}
TESORTER_ENV=${TESORTER_ENV:-tesorter}

activate_edta() {
    local status=0 had_nounset=0
    case $- in *u*) had_nounset=1;; esac
    # Conda compiler hooks may read unset backup variables during an environment switch.
    set +u
    source "$CONDA_SH" && conda activate "$EDTA_ENV" || status=$?
    if [ "$had_nounset" -eq 1 ]; then set -u; fi
    return "$status"
}

activate_tesorter() {
    local status=0 had_nounset=0
    case $- in *u*) had_nounset=1;; esac
    # Conda compiler hooks may read unset backup variables during an environment switch.
    set +u
    source "$CONDA_SH" && conda activate "$TESORTER_ENV" || status=$?
    if [ "$had_nounset" -eq 1 ]; then set -u; fi
    return "$status"
}

count_matches() {
    local count status
    if count=$(grep -cE -- "$1" "$2"); then
        :
    else
        status=$?
        [ "$status" -eq 1 ] || return "$status"
    fi
    printf '%s\n' "$count"
}

require_success() {
    [ "$2" -eq 0 ] && return 0
    printf 'ERROR: %s exited with code %s\n' "$1" "$2" >&2
    return "$2"
}

select_sample() {
    local task_id=${SLURM_ARRAY_TASK_ID:-1}
    case "$task_id" in ''|*[!0-9]*) echo 'Invalid task ID' >&2; return 1;; esac
    [ "$task_id" -ge 1 ] || return 1
    SAMPLE=$(awk -v n="$task_id" 'NR==n {print; exit}' "$SAMPLE_TABLE") || return 1
    [ -n "$SAMPLE" ] || { echo 'Sample task ID out of range' >&2; return 1; }
}

select_chromosome() {
    local task_id=${SLURM_ARRAY_TASK_ID:-1} row
    case "$task_id" in ''|*[!0-9]*) echo 'Invalid task ID' >&2; return 1;; esac
    [ "$task_id" -ge 1 ] || return 1
    row=$(awk -v n="$task_id" 'NR==n {print; exit}' "$TASK_TABLE") || return 1
    IFS=$'\t' read -r SAMPLE CHR <<< "$row"
    [ -n "$SAMPLE" ] && [ -n "$CHR" ] || { echo 'Chromosome task ID out of range' >&2; return 1; }
}

edta_cds_source() {
    local cds="$CDS_DIR/$1.cds.fa"
    if [ -f "$cds" ]; then printf '%s\n' "$cds"; else printf 'none\n'; fi
}

record_edta_success() {
    local marker="$2/edta.complete"
    printf 'complete\n%s\n' "$(edta_cds_source "$1")" > "$marker.tmp" && mv "$marker.tmp" "$marker"
}

edta_log_clean() {
    local failures
    failures=$(count_matches "daemonic processes are not allowed|TIR results not found|SINE results not found|Helitron.*(没有那个文件|No such file)|Can't locate Bio|LTR_FINDER_parallel.*No such" "$1") || return 1
    [ "$failures" -eq 0 ]
}

current_chromosome_result() {
    local sample=$1 chr=$2 genome work library log marker cds input
    genome="$GENOME_DIR/$sample.chrom.fasta"
    work="$EDTA_OUT_DIR/$sample/per_chr/$chr"
    library="$work/$chr.fa.mod.EDTA.TElib.fa"
    log="$work/edta_$chr.log"
    marker="$work/edta.complete"
    for input in "$genome" "$genome.fai" "$library" "$log" "$marker"; do
        [ -s "$input" ] || return 1
        [ "$input" -nt "$marker" ] && return 1
    done
    cds=$(edta_cds_source "$sample") || return 1
    [ "$(sed -n '2p' "$marker")" = "$cds" ] || return 1
    [ "$cds" != none ] && [ "$cds" -nt "$marker" ] && return 1
    edta_log_clean "$log"
}

current_reclassified_library() {
    local sample=$1 genome fai merged nr library nr_mark lib_mark chr rest chr_dir chr_lib chr_mark
    genome="$GENOME_DIR/$sample.chrom.fasta"
    fai="$genome.fai"
    merged="$EDTA_OUT_DIR/$sample/merged"
    nr="$merged/${sample}_TElib_nr.fa"
    library="$merged/${sample}_TElib_nr.reclassified.fa"
    nr_mark="$merged/nr.complete"
    lib_mark="$merged/reclass.complete"
    [ -s "$genome" ] && [ -s "$fai" ] && [ -s "$nr" ] && [ -s "$library" ] && [ -f "$nr_mark" ] && [ -f "$lib_mark" ] || return 1
    [ "$genome" -nt "$nr_mark" ] || [ "$fai" -nt "$nr_mark" ] || [ "$nr" -nt "$nr_mark" ] || [ "$library" -nt "$lib_mark" ] || [ "$nr" -nt "$lib_mark" ] || [ "$nr_mark" -nt "$lib_mark" ] && return 1
    while read -r chr rest; do
        chr_dir="$EDTA_OUT_DIR/$sample/per_chr/$chr"
        chr_lib="$chr_dir/$chr.fa.mod.EDTA.TElib.fa"
        chr_mark="$chr_dir/edta.complete"
        current_chromosome_result "$sample" "$chr" || return 1
        [ "$chr_mark" -nt "$nr_mark" ] || [ "$chr_lib" -nt "$nr_mark" ] && return 1
    done < "$fai"
    return 0
}
