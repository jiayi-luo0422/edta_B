#!/bin/bash
#SBATCH --job-name=rm_panTE
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=normal
#SBATCH --mem=100G
#SBATCH --time=120:00:00
#SBATCH --array=1
#SBATCH --output=logs/rm_panTE_%A_%a.out
#SBATCH --error=logs/rm_panTE_%A_%a.err

# 06: 用第05步最终pan-TE库统一重注释样本表中的全部组装。
# 数组范围在提交时按SAMPLE_TABLE行数设置；结果独立于04自库注释。

set -euo pipefail
SCRIPT_DIR=$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)
source "${SCRIPT_DIR}/workflow_common.sh"
activate_edta

THREADS=4
select_sample

GENOME="${GENOME_DIR}/${SAMPLE}.chrom.fasta"
TELIB="${PANTE_DIR}/Tea_panTE.lib.fasta"
RM_DIR="${PANTE_MASK_DIR}/${SAMPLE}"
RM_OUT="${RM_DIR}/$(basename "${GENOME}").out"
RM_MARK="${RM_DIR}/repeatmasker.complete"

echo "============================================================"
echo "[$(date)] RepeatMasker (final pan-TE lib): ${SAMPLE}"
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
if [ ! -s "${TELIB}" ] || ! grep -q '^>' "${TELIB}"; then
    echo "  ✗ 最终pan-TE库缺失或无FASTA记录: ${TELIB}"
    echo "    请先完成05_build_pan_te_library.sh并验收"
    exit 1
fi

# 记录输入来源，库被替换或切换路径时不能复用旧注释；不计算文件哈希。
input_record() {
    python3 - "${GENOME}" "${TELIB}" <<'PYRECORD'
import json
from pathlib import Path
import sys

records = []
for name in sys.argv[1:]:
    path = Path(name).resolve()
    stat = path.stat()
    records.append({'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns})
print(json.dumps({'inputs': records, 'parameters': ['-pa', '4', '-gff', '-no_is', '-nolow', '-xsmall']}, sort_keys=True))
PYRECORD
}
INPUT_RECORD=$(input_record)

########################################################################
# 3. 断点续跑：成功标记、完整产物与输入记录
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
    # 逐行流式比较，忽略序列大小写与换行；名称、顺序和碱基必须一致。
    cmp -s <(normalized_fasta "${GENOME}") <(normalized_fasta "${base}.masked")
}

normalized_fasta() {
    awk '/^>/{if(NR>1)printf "\n"; print $1; next}
         {gsub(/[[:space:]]/, ""); printf "%s", toupper($0)}
         END{printf "\n"}' "$1"
}

if [ -f "${RM_MARK}" ] && [ -f "${RM_DIR}/repeatmasker.inputs.json" ] &&
   cmp -s <(printf '%s\n' "${INPUT_RECORD}") "${RM_DIR}/repeatmasker.inputs.json" &&
   validate_rm_outputs "${RM_DIR}"; then
    echo "  [跳过] ${SAMPLE} 已验收且输入未更新"
    exit 0
fi

mkdir -p "$(dirname "${RM_DIR}")"
RUN_DIR=$(mktemp -d "${RM_DIR}.run.XXXXXX")
printf '%s\n' "${INPUT_RECORD}" > "${RUN_DIR}/repeatmasker.inputs.json"

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
if [ "${INPUT_RECORD}" != "$(input_record)" ]; then
    echo "  ✗ 运行期间基因组或pan-TE库改变；旧结果保留，检查 ${RUN_DIR}"; exit 1
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
echo "  → .out 使用统一pan-TE库中的标签；本流程完成"
