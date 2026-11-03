#!/bin/bash
#SBATCH --job-name=edta_v2
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=30:00:00
#SBATCH --array=1
#SBATCH --output=logs/edta_v2_%A_%a.out
#SBATCH --error=logs/edta_v2_%A_%a.err

#==========================================================================
# 02_annotate_chromosomes.sh
# 按 edta_tasks.tsv 每行的 sample / chromosome 提交一个任务。
# 先运行 01_prepare_task_table.sh；数组范围取该表行数，详见 README.md。
# 已完成且 edta.complete 与输入保持一致的染色体可复用。
# 合并按 edta_samples.txt 的行数提交 03_merge_and_classify_libraries.sh。
#==========================================================================

set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/workflow_common.sh" || exit 1
activate_edta || exit 1

# ----- 配置 -----
THREADS=16

# Explicit task table prepared by 01_prepare_task_table.sh.
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
select_chromosome || exit 1
GENOME_ORIG="${GENOME_DIR}/${SAMPLE}.chrom.fasta"
GENOME_FAI="${GENOME_ORIG}.fai"
[ -s "$GENOME_ORIG" ] && [ -s "$GENOME_FAI" ] || { echo "Missing genome/index: $GENOME_ORIG" >&2; exit 1; }

WORKDIR="${EDTA_OUT_DIR}/${SAMPLE}/per_chr/${CHR}"
CDS_FILE="${CDS_DIR}/${SAMPLE}.cds.fa"

echo "============================================================"
echo "[$(date)] EDTA v2 注释: task ${TASK_ID} → ${SAMPLE}/${CHR}"
echo "  WORKDIR: ${WORKDIR}"
echo "  env:     EDTA2_fix"
echo "  CPU/MEM: ${THREADS} / 100G"
echo "============================================================"

########################################################################
# 1. 断点续跑检查
########################################################################
TELIB_EXPECT="${WORKDIR}/${CHR}.fa.mod.EDTA.TElib.fa"
LOG_EXPECT="${WORKDIR}/edta_${CHR}.log"

if [ -s "${TELIB_EXPECT}" ]; then
    if ! current_chromosome_result "$SAMPLE" "$CHR"; then
        echo "Existing library has no current success record; inspect it before using optional/rerun_chromosome.sh: $WORKDIR" >&2
        exit 1
    fi
    echo "  [跳过] 已完成: $SAMPLE/$CHR"
    exit 0
fi

mkdir -p "${WORKDIR}" || exit 1
cd "${WORKDIR}" || exit 1

########################################################################
# 2. 提取染色体序列
########################################################################
CHR_FASTA="${WORKDIR}/${CHR}.fa"
if [ ! -s "${CHR_FASTA}" ] || [ "$GENOME_ORIG" -nt "$CHR_FASTA" ]; then
    echo "[$(date)] 提取 ${CHR}"
    "$SAMTOOLS_BIN" faidx "${GENOME_ORIG}" "${CHR}" > "${CHR_FASTA}" || exit 1
    [ -s "${CHR_FASTA}" ] || exit 1
    SIZE_MB=$(awk -v s="$(stat --printf='%s' "${CHR_FASTA}")" 'BEGIN{printf "%.0f", s/1024/1024}')
    echo "  大小: ${SIZE_MB} MB"
fi

########################################################################
# 3. 跑 EDTA
########################################################################
echo ""
echo "[$(date)] 启动 EDTA"

EDTA_CMD=(EDTA.pl \
    --genome "${CHR_FASTA}" \
    --species others \
    --step all \
    --anno 0 \
    --sensitive 0 \
    --force 1 \
    --threads "${THREADS}")

if [ -f "${CDS_FILE}" ]; then
    EDTA_CMD+=(--cds "${CDS_FILE}")
fi

rm -f "${WORKDIR}/edta.complete" || exit 1
START_EPOCH=$(date +%s)
"${EDTA_CMD[@]}" 2>&1 | tee "${LOG_EXPECT}"
EDTA_STATUS=("${PIPESTATUS[@]}")
EDTA_EXIT=${EDTA_STATUS[0]}
require_success "EDTA log writer" "${EDTA_STATUS[1]}" || exit 1
END_EPOCH=$(date +%s)
ELAPSED_HOURS=$(awk -v s="${START_EPOCH}" -v e="${END_EPOCH}" 'BEGIN{printf "%.1f", (e-s)/3600}')

echo ""
echo "[$(date)] ${SAMPLE}/${CHR} EDTA 退出码: ${EDTA_EXIT}, 耗时 ${ELAPSED_HOURS}h"

########################################################################
# 4. 自检
########################################################################
n_daemonic=$(count_matches "daemonic processes are not allowed" "${LOG_EXPECT}") || exit 1
n_tir_missing=$(count_matches "TIR results not found" "${LOG_EXPECT}") || exit 1

if [ -s "${TELIB_EXPECT}" ]; then
    n_te=$(grep -c "^>" "${TELIB_EXPECT}")
else
    n_te="MISSING"
fi

echo ""
echo "  自检:"
echo "    daemonic 错误:    ${n_daemonic}"
echo "    TIR missing:      ${n_tir_missing}"
echo "    TElib families:   ${n_te}"

if ! edta_log_clean "$LOG_EXPECT"; then
    echo "  ✗ ${SAMPLE}/${CHR} 仍有失败信号"
    exit 2
fi

if [ "${EDTA_EXIT}" -ne 0 ] || [ ! -s "${TELIB_EXPECT}" ]; then
    echo "  ✗ ${SAMPLE}/${CHR} EDTA 未正常完成"
    exit 3
fi

record_edta_success "$SAMPLE" "$WORKDIR" || exit 1
echo "  ✓ ${SAMPLE}/${CHR} 任务通过"

########################################################################
# 5. 清理大型中间产物 (释放磁盘)
#    保留: TElib.fa, intact.gff3, log
#    删除: raw 目录下大量 BLAST 比对产物
########################################################################
echo ""
echo "[$(date)] 清理中间产物以省磁盘 (保留 TElib + intact + log)"

# 这些是大文件且后续不需要
RAW_DIR="${CHR_FASTA}.mod.EDTA.raw"
ANNO_DIR="${CHR_FASTA}.mod.EDTA.anno"
COMBINE_DIR="${CHR_FASTA}.mod.EDTA.combine"
FINAL_DIR="${CHR_FASTA}.mod.EDTA.final"

# 保留 raw 目录里的 *.intact.raw.fa 和 *.raw.fa, 删 BLAST 中间数据
if [ -d "${RAW_DIR}" ]; then
    find "${RAW_DIR}" -type f \
        ! -name "*.intact.raw.fa" \
        ! -name "*.raw.fa" \
        ! -name "*.gff3" \
        ! -name "*.bed" \
        ! -name "*.list" \
        -delete 2>/dev/null || true
    # 删空子目录
    find "${RAW_DIR}" -type d -empty -delete 2>/dev/null || true
fi

# 这两个目录后续完全不需要
rm -rf "${COMBINE_DIR}" 2>/dev/null || true

# .mod 文件本身和 RM2.fa, MAKER.masked 等也大, 但保留 .mod.EDTA.TElib.fa (核心输出)
# 删除巨大但下游不要的:
rm -f "${CHR_FASTA}.mod.MAKER.masked" 2>/dev/null || true
rm -f "${CHR_FASTA}.mod.EDTA.TEanno.gff3" 2>/dev/null || true   # --anno 0 时也可能产生
rm -f "${CHR_FASTA}.mod.EDTA.TEanno.density_plots.pdf" 2>/dev/null || true

# 显示清理后的目录大小
SIZE_AFTER=$(du -sh "${WORKDIR}" 2>/dev/null | cut -f1)
echo "  ${WORKDIR} 清理后大小: ${SIZE_AFTER}"

echo ""
echo "============================================================"
echo "[$(date)] task ${TASK_ID} (${SAMPLE}/${CHR}) 完成"
echo "============================================================"
