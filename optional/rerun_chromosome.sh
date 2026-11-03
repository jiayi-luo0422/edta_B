#!/bin/bash
#SBATCH --job-name=edta_rescue
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --output=logs/rescue_%j.out
#SBATCH --error=logs/rescue_%j.err

#==========================================================================
# optional/rerun_chromosome.sh
# 通用单染色体 EDTA 救援/补丁脚本
#
# 使用方法:
#   1. 设置 RERUN_SAMPLE 和 RERUN_CHR 为需重跑的样本和染色体
#   2. 在根脚本目录投递: sbatch optional/rerun_chromosome.sh
#==========================================================================

set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "$(dirname "${BASH_SOURCE[0]}")")}/workflow_common.sh" || exit 1
activate_edta || exit 1

## =================== 参数配置区 ===================
SAMPLE="${RERUN_SAMPLE:?Set RERUN_SAMPLE}"
CHR="${RERUN_CHR:?Set RERUN_CHR}"
## ===================================================

THREADS=8

WORKDIR="${EDTA_OUT_DIR}/${SAMPLE}/per_chr/${CHR}"
CDS_FILE="${CDS_DIR}/${SAMPLE}.cds.fa"
GENOME_ORIG="${GENOME_DIR}/${SAMPLE}.chrom.fasta"
CHR_FASTA="${WORKDIR}/${CHR}.fa"
LOG_EXPECT="${WORKDIR}/edta_${CHR}.log"
TELIB_EXPECT="${WORKDIR}/${CHR}.fa.mod.EDTA.TElib.fa"

mkdir -p "${WORKDIR}" || exit 1
cd "${WORKDIR}" || exit 1

echo "============================================================"
echo "[$(date)] EDTA 独立补丁启动: ${SAMPLE} / ${CHR}"
echo "  WORKDIR: ${WORKDIR}"
echo "  资源配置: ${THREADS} CPU / 150G Mem"
echo "============================================================"

########################################################################
# 1. 提取序列检查
########################################################################
if [ ! -s "${CHR_FASTA}" ] || [ "$GENOME_ORIG" -nt "$CHR_FASTA" ]; then
    echo "[$(date)] 提取染色体序列: ${CHR}"
    "$SAMTOOLS_BIN" faidx "${GENOME_ORIG}" "${CHR}" > "${CHR_FASTA}" || exit 1
    [ -s "${CHR_FASTA}" ] || exit 1
fi

########################################################################
# 2. 核心步骤：强力清理旧的损坏缓存（防止 Silent Failure 和死锁残留）
########################################################################
RAW_DIR="${CHR_FASTA}.mod.EDTA.raw"
if [ -d "${RAW_DIR}" ]; then
    echo "[$(date)] 检测到历史残缺目录，正在强制清理以释放缓存: ${RAW_DIR}"
    rm -rf "${RAW_DIR}"
fi

# 顺便清理可能残留的最终生成标志，确保本次跑出来的是纯净的新结果
rm -f "${TELIB_EXPECT}"

########################################################################
# 3. 运行 EDTA
########################################################################
echo ""
echo "[$(date)] 启动 EDTA.pl 主程序..."

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
echo "[$(date)] EDTA 运行完毕，退出码: ${EDTA_EXIT}，耗时: ${ELAPSED_HOURS} 小时"

########################################################################
# 4. 自动化自检
########################################################################
echo ""
echo "[$(date)] 开始运行结果自检..."
n_daemonic=$(count_matches "daemonic processes are not allowed" "${LOG_EXPECT}") || exit 1
n_sine_missing=$(count_matches "SINE results not found" "${LOG_EXPECT}") || exit 1
n_tir_missing=$(count_matches "TIR results not found" "${LOG_EXPECT}") || exit 1

if [ -s "${TELIB_EXPECT}" ]; then
    n_te=$(grep -c "^>" "${TELIB_EXPECT}")
else
    n_te="MISSING"
fi

echo "  自检报告:"
echo "    - TElib 产生条数:  ${n_te}"
echo "    - Deadlock 错误数: ${n_daemonic}"
echo "    - SINE 缺失信号:   ${n_sine_missing}"
echo "    - TIR 缺失信号:    ${n_tir_missing}"

if ! edta_log_clean "$LOG_EXPECT" || [ "${EDTA_EXIT}" -ne 0 ] || [ "${n_te}" = "MISSING" ] || [ "${n_daemonic}" -gt 0 ] || [ "${n_sine_missing}" -gt 0 ] || [ "${n_tir_missing}" -gt 0 ]; then
    echo "  ✗ 自检失败！该染色体仍存在异常，请检查 log。"
    exit 1
else
    record_edta_success "$SAMPLE" "$WORKDIR" || exit 1
    echo "  ✓ 计算与产物检查通过"
fi

########################################################################
# 5. 按照集群标准清理大型中间产物（释放磁盘空间）
########################################################################
echo ""
echo "[$(date)] 清理低价值中间比对数据..."

COMBINE_DIR="${CHR_FASTA}.mod.EDTA.combine"
if [ -d "${RAW_DIR}" ]; then
    find "${RAW_DIR}" -type f \
        ! -name "*.intact.raw.fa" \
        ! -name "*.raw.fa" \
        ! -name "*.gff3" \
        ! -name "*.bed" \
        ! -name "*.list" \
        -delete 2>/dev/null || true
    find "${RAW_DIR}" -type d -empty -delete 2>/dev/null || true
fi

rm -rf "${COMBINE_DIR}" 2>/dev/null || true
rm -f "${CHR_FASTA}.mod.MAKER.masked" 2>/dev/null || true
rm -f "${CHR_FASTA}.mod.EDTA.TEanno.gff3" 2>/dev/null || true
rm -f "${CHR_FASTA}.mod.EDTA.TEanno.density_plots.pdf" 2>/dev/null || true

SIZE_AFTER=$(du -sh "${WORKDIR}" 2>/dev/null | cut -f1)
echo "  工作目录清理后最终大小: ${SIZE_AFTER}"
echo "============================================================"