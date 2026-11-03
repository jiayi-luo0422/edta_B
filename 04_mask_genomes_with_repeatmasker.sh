#!/bin/bash
#SBATCH --job-name=rm_per_sample
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=normal
#SBATCH --mem=100G
#SBATCH --time=120:00:00
#SBATCH --array=1
#SBATCH --output=logs/rm_per_sample_%A_%a.out
#SBATCH --error=logs/rm_per_sample_%A_%a.err

#==========================================================================
# 04_mask_genomes_with_repeatmasker.sh   (干净版)
# 每品种用自家重分类库跑 RepeatMasker -> .out
#
#
# 前置: 03 (含 TEsorter 重分类) 完成
# 完成后: sbatch 05_build_pan_te_library.sh
#==========================================================================

set -euo pipefail
SCRIPT_DIR=$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)
source "${SCRIPT_DIR}/workflow_common.sh"
activate_edta

THREADS=4
select_sample

GENOME="${GENOME_DIR}/${SAMPLE}.chrom.fasta"
TELIB="${EDTA_OUT_DIR}/${SAMPLE}/merged/${SAMPLE}_TElib_nr.reclassified.fa"   # ★ 只认重分类库
RM_DIR="${EDTA_OUT_DIR}/${SAMPLE}/repeatmasker_self"
RM_OUT="${RM_DIR}/$(basename "${GENOME}").out"
RM_MARK="${RM_DIR}/repeatmasker.complete"

echo "============================================================"
echo "[$(date)] RepeatMasker (self reclassified lib): ${SAMPLE}"
echo "  genome: ${GENOME}"
echo "  lib:    ${TELIB}"
echo "  out:    ${RM_DIR}"
echo "============================================================"

########################################################################
# 2. 前置检查
########################################################################
if [ ! -s "${GENOME}" ]; then
    echo "  ✗ genome 缺失: ${GENOME}"; exit 1
fi
if [ ! -s "${TELIB}" ]; then
    echo "  ✗ ${SAMPLE} 重分类库缺失: ${TELIB}"
    echo "    请先跑 03_merge_and_classify_libraries.sh (含 TEsorter 重分类)"
    exit 1
fi

if ! current_reclassified_library "${SAMPLE}"; then
    echo "  ✗ ${SAMPLE} 重分类库未验收或上游输入已更新，请先完成02/03"; exit 1
fi

########################################################################
# 3. 断点续跑：成功标记、完整产物与两个输入的时间戳
########################################################################
zero_rm_output() {
    awk 'NR==1 && /^There were no repetitive sequences detected in .+/ {zero=1}
         END {exit !(zero && NR==1)}' "$1"
}

validate_rm_outputs() {
    local base="$1/$(basename "${GENOME}")"
    [ -s "${base}.out" ] && [ -s "${base}.masked" ] || return 1
    # ProcessRepeats 的零命中分支不生成常规表格或 GFF。
    if ! zero_rm_output "${base}.out"; then
        [ -s "${base}.tbl" ] && [ -f "${base}.out.gff" ] || return 1
        awk 'NR<=3 && /SW[[:space:]]+perc/ {header=1} END {exit !header}' "${base}.out" || return 1
    fi
    # Masked FASTA 必须保留基因组的序列名和每条序列长度；允许合法零TE结果。
    cmp -s <(awk '/^>/{if(n!="")print n,len; n=$1; len=0; next} {gsub(/[[:space:]]/, ""); len+=length} END{if(n!="")print n,len}' "${GENOME}") \
           <(awk '/^>/{if(n!="")print n,len; n=$1; len=0; next} {gsub(/[[:space:]]/, ""); len+=length} END{if(n!="")print n,len}' "${base}.masked")
}

if [ -f "${RM_MARK}" ] && [ ! "${GENOME}" -nt "${RM_MARK}" ] &&
   [ ! "${TELIB}" -nt "${RM_MARK}" ] && [ ! "${GENOME}" -nt "${RM_OUT}" ] &&
   [ ! "${TELIB}" -nt "${RM_OUT}" ] && validate_rm_outputs "${RM_DIR}"; then
    echo "  [跳过] ${SAMPLE} 已验收且输入未更新"
    exit 0
fi

mkdir -p "$(dirname "${RM_DIR}")"
RUN_DIR=$(mktemp -d "${RM_DIR}.run.XXXXXX")

########################################################################
# 4. 跑 RepeatMasker
########################################################################
START_EPOCH=$(date +%s)
echo ""
echo "[$(date)] 启动 RepeatMasker"

set +e
RepeatMasker \
    -pa ${THREADS} \
    -lib "${TELIB}" \
    -dir "${RUN_DIR}" \
    -gff -no_is -nolow -xsmall \
    "${GENOME}" \
    2>&1 | tee "${RUN_DIR}/repeatmasker.log"
RM_STATUS=("${PIPESTATUS[@]}")
set -e
RM_EXIT=${RM_STATUS[0]}
[ "${RM_STATUS[1]}" -eq 0 ] || RM_EXIT=${RM_STATUS[1]}
END_EPOCH=$(date +%s)
ELAPSED_HOURS=$(awk -v s="${START_EPOCH}" -v e="${END_EPOCH}" 'BEGIN{printf "%.1f", (e-s)/3600}')
echo ""
echo "[$(date)] ${SAMPLE} RepeatMasker 退出码: ${RM_EXIT}, 耗时 ${ELAPSED_HOURS}h"

########################################################################
# 5. 自检
########################################################################
require_success "RepeatMasker (${RUN_DIR}/repeatmasker.log；旧结果保留)" "${RM_EXIT}"
# 零命中且退出成功时，原始 FASTA 即为等价的软屏蔽结果。
if zero_rm_output "${RUN_DIR}/$(basename "${GENOME}").out" &&
   [ ! -e "${RUN_DIR}/$(basename "${GENOME}").masked" ]; then
    cp "${GENOME}" "${RUN_DIR}/$(basename "${GENOME}").masked"
fi
if ! validate_rm_outputs "${RUN_DIR}"; then
    echo "  ✗ RepeatMasker产物不完整或masked序列不完整；旧结果保留，检查 ${RUN_DIR}"; exit 1
fi
# 所有检查通过后才发布；旧目录在失败的重跑中保持可用。
printf "complete\n" > "${RUN_DIR}/repeatmasker.complete"
if [ -d "${RM_DIR}" ]; then
    OLD_DIR="${RUN_DIR}.previous"
    mv "${RM_DIR}" "${OLD_DIR}"
    if ! mv "${RUN_DIR}" "${RM_DIR}"; then
        mv "${OLD_DIR}" "${RM_DIR}"; exit 1
    fi
    rm -rf "${OLD_DIR}"
else
    mv "${RUN_DIR}" "${RM_DIR}"
fi
n_lines=$(wc -l < "${RM_OUT}")
n_te=$(awk 'NR>3 && NF>10' "${RM_OUT}" | wc -l)
echo ""
echo "  ✓ ${SAMPLE}.out: ${n_lines} 行 (${n_te} TE 注释)"
echo "  → .out 使用重分类库中的标签；下一步见 05_build_pan_te_library.sh"
