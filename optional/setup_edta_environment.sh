#!/bin/bash
#SBATCH --job-name=edta_setup
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/edta_setup_%j.out
#SBATCH --error=logs/edta_setup_%j.err

#==========================================================================
# optional/setup_edta_environment.sh
#
# 创建 Python 3.10 + EDTA 2.2.2 的新 conda env, 同时修复:
#   1. 配置指定 Python / EDTA 版本
#   2. BioPerl 缺失 (LTR_retriever 静默失败)         ← 历史坑 1
#   3. LTR_FINDER_parallel / LTR_HARVEST_parallel 路径错位  ← 历史坑 2
#
# 在根脚本目录提交安装作业:
#   sbatch optional/setup_edta_environment.sh
#
# 完成后:
#   - 新 env: EDTA2_fix
#   - 验证文件: ~/EDTA2_fix_validation.txt
#   - 源码备份: ~/EDTA2_fix_pristine_backup/
#==========================================================================

set -uo pipefail

source "${SLURM_SUBMIT_DIR:-$(dirname "$(dirname "${BASH_SOURCE[0]}")")}/workflow_common.sh" || exit 1
ENV_NAME="$EDTA_ENV"
PYTHON_VER="3.10"
EDTA_VER="2.2.2"
VALIDATION="${HOME}/EDTA2_fix_validation.txt"

echo "============================================================"
echo "EDTA2_fix 环境创建 v2  $(date)"
echo "============================================================"

########################################################################
# 1. 加载 conda
########################################################################
source "$CONDA_SH" || exit 1

########################################################################
# 2. 检查并创建环境
########################################################################
if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "  ⚠ env ${ENV_NAME} 已存在"
    echo "  如需重建, 先 conda env remove -n ${ENV_NAME} -y, 然后再跑此脚本"
    echo "  当前直接进入修复 + 验证流程"
else
    echo "  [$(date)] 创建 ${ENV_NAME}: python=${PYTHON_VER} edta=${EDTA_VER} + BioPerl"
    echo "  (使用 mamba 加速, 同时显式装 perl-bioperl 避开历史坑 1)"

    # 先试 mamba, 失败回退 conda
    if command -v mamba >/dev/null 2>&1; then
        SOLVER="mamba"
    else
        SOLVER="conda"
        echo "  注意: mamba 不可用, 用 conda 解依赖 (慢一些)"
    fi

    ${SOLVER} create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        python="${PYTHON_VER}" \
        edta="${EDTA_VER}" \
        perl-bioperl \
        samtools \
        cd-hit \
        mmseqs2 \
        repeatmasker \
        -y 2>&1 | tail -30

    if [ $? -ne 0 ]; then
        echo "  ✗ 环境创建失败"
        exit 1
    fi
fi

########################################################################
# 3. 激活
########################################################################
echo ""
echo "[$(date)] 激活 ${ENV_NAME}"
activate_edta || exit 1

# 锚点: EDTA 安装根 (这一行很重要, 后续所有路径修复都基于这里)
EDTA_PREFIX=$(dirname $(which EDTA.pl))/..
EDTA_PREFIX=$(realpath "${EDTA_PREFIX}")
EDTA_BASE="${EDTA_PREFIX}/share/EDTA"

echo "  EDTA_PREFIX: ${EDTA_PREFIX}"
echo "  EDTA_BASE:   ${EDTA_BASE}"

########################################################################
# 4. 修复历史坑 2: LTR_FINDER_parallel / LTR_HARVEST_parallel 路径错位
########################################################################
echo ""
echo "[$(date)] 修复 LTR_FINDER/HARVEST_parallel 路径软链"

# 4a. LTR_FINDER_parallel
if [ -f "${EDTA_BASE}/bin/LTR_FINDER_parallel/LTR_FINDER_parallel" ]; then
    echo "  ✓ LTR_FINDER_parallel 已就位"
else
    SRC1="${EDTA_PREFIX}/share/LTR_FINDER_parallel/LTR_FINDER_parallel"
    if [ -f "${SRC1}" ]; then
        mkdir -p "${EDTA_BASE}/bin/LTR_FINDER_parallel"
        ln -sf "${SRC1}" "${EDTA_BASE}/bin/LTR_FINDER_parallel/LTR_FINDER_parallel"
        echo "  ✓ 建立 LTR_FINDER_parallel 软链"
    else
        echo "  ⚠ 找不到 LTR_FINDER_parallel 源文件: ${SRC1}"
        echo "    手动查找:"
        find "${EDTA_PREFIX}/share" -name "LTR_FINDER_parallel" -type f 2>/dev/null | head -3
    fi
fi

# FindBin uses the entry directory, so its helper directory must also be reachable.
FINDER_DIR="${EDTA_BASE}/bin/LTR_FINDER_parallel"
if [ ! -e "${FINDER_DIR}/bin" ] && [ ! -L "${FINDER_DIR}/bin" ]; then
    FINDER_BIN="${EDTA_PREFIX}/share/LTR_FINDER_parallel/bin"
    [ -d "${FINDER_BIN}" ] || { echo "Missing LTR_FINDER helpers: ${FINDER_BIN}" >&2; exit 1; }
    ln -s "${FINDER_BIN}" "${FINDER_DIR}/bin" || exit 1
fi
for helper in cut.pl convert_ltr_finder2.pl; do
    [ -r "${FINDER_DIR}/bin/${helper}" ] && [ -s "${FINDER_DIR}/bin/${helper}" ] || {
        echo "Missing LTR_FINDER helper: ${FINDER_DIR}/bin/${helper}" >&2
        exit 1
    }
done

# 4b. LTR_HARVEST_parallel
if [ -f "${EDTA_BASE}/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel" ]; then
    echo "  ✓ LTR_HARVEST_parallel 已就位"
else
    SRC2="${EDTA_PREFIX}/share/LTR_HARVEST_parallel/LTR_HARVEST_parallel"
    if [ -f "${SRC2}" ]; then
        mkdir -p "${EDTA_BASE}/bin/LTR_HARVEST_parallel"
        ln -sf "${SRC2}" "${EDTA_BASE}/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel"
        echo "  ✓ 建立 LTR_HARVEST_parallel 软链"
    else
        echo "  ⚠ 找不到 LTR_HARVEST_parallel 源文件: ${SRC2}"
        find "${EDTA_PREFIX}/share" -name "LTR_HARVEST_parallel" -type f 2>/dev/null | head -3
    fi
fi

########################################################################
# 5. 全量验证 (写入 ~/EDTA2_fix_validation.txt)
########################################################################
echo ""
echo "[$(date)] 全量验证"

> "${VALIDATION}"
{
    echo "EDTA2_fix 环境验证报告 v2"
    echo "生成时间: $(date)"
    echo "============================================================"
    echo ""

    echo "[A] Python 与 EDTA 版本"
    python --version 2>&1
    which python
    EDTA.pl --help 2>&1 | head -3
    which EDTA.pl
    echo "EDTA_PREFIX: ${EDTA_PREFIX}"
    echo ""

    echo "[B] Python 版本与现有模块（TIR-Learner 4不要求dask/swifter）"
    python -c "import sys; print(f'Python: {sys.version}')"
    python -c "import pandas; print(f'pandas: {pandas.__version__}')" 2>&1
    python -c "import Bio; print(f'biopython: {Bio.__version__}')" 2>&1
    echo ""

    echo "[C] BioPerl (历史坑 1)"
    perl -e "use Bio::SeqIO; print 'Bio::SeqIO loaded\n';" 2>&1
    echo ""

    echo "[D] LTR_FINDER/HARVEST_parallel 路径 (历史坑 2)"
    if [ -e "${EDTA_BASE}/bin/LTR_FINDER_parallel/LTR_FINDER_parallel" ]; then
        echo "  ✓ LTR_FINDER_parallel: $(readlink -f ${EDTA_BASE}/bin/LTR_FINDER_parallel/LTR_FINDER_parallel)"
    else
        echo "  ✗ LTR_FINDER_parallel 缺失"
    fi
    if [ -e "${EDTA_BASE}/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel" ]; then
        echo "  ✓ LTR_HARVEST_parallel: $(readlink -f ${EDTA_BASE}/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel)"
    else
        echo "  ✗ LTR_HARVEST_parallel 缺失"
    fi
    echo ""

    echo "[E] EDTA 内部脚本路径全扫描"
    if [ -f "${EDTA_BASE}/EDTA_raw.pl" ]; then
        for script in $(grep -ohP 'bin/[^/"]+/[^/"]+' \
            $(which EDTA.pl) "${EDTA_BASE}/EDTA_raw.pl" 2>/dev/null | sort -u); do
            path="${EDTA_BASE}/${script}"
            if [ -e "$path" ]; then
                echo "  OK     : $script"
            else
                echo "  缺失   : $script"
            fi
        done
    else
        echo "  跳过(EDTA_raw.pl 不存在)"
    fi
    echo ""

    echo "[F] 核心工具命令"
    for tool in ltr_finder gt TIR-Learner HelitronScanner RepeatModeler RepeatMasker samtools cd-hit-est mmseqs; do
        path=$(which "${tool}" 2>/dev/null)
        if [ -n "${path}" ]; then
            echo "  ✓ ${tool}: ${path}"
        else
            echo "  ✗ ${tool}: 未找到"
        fi
    done
    echo ""


} | tee "${VALIDATION}"

########################################################################
# 6. 备份关键文件 (为后续可能的补丁做准备)
########################################################################
echo ""
echo "[$(date)] 备份 EDTA 关键文件"

BACKUP_DIR="${HOME}/EDTA2_fix_pristine_backup"
mkdir -p "${BACKUP_DIR}"

if [ ! -d "${BACKUP_DIR}/EDTA" ] && [ -d "${EDTA_BASE}" ]; then
    cp -rL "${EDTA_BASE}" "${BACKUP_DIR}/" 2>/dev/null || \
    cp -r  "${EDTA_BASE}" "${BACKUP_DIR}/" 2>/dev/null
    echo "  ✓ EDTA 源码备份 → ${BACKUP_DIR}/EDTA"
fi

TIRLEARNER_DIR=$(find "${EDTA_PREFIX}/share" -maxdepth 2 -name "TIR-Learner*" -type d 2>/dev/null | head -1)
if [ -n "${TIRLEARNER_DIR}" ] && [ ! -d "${BACKUP_DIR}/$(basename ${TIRLEARNER_DIR})" ]; then
    cp -r "${TIRLEARNER_DIR}" "${BACKUP_DIR}/"
    echo "  ✓ TIR-Learner 源码备份 → ${BACKUP_DIR}/$(basename ${TIRLEARNER_DIR})"
fi

########################################################################
# 7. 最终总结
########################################################################
echo ""
echo "============================================================"
echo "环境创建完成"
echo "============================================================"
echo "  env:           ${ENV_NAME}"
echo "  Python:        ${PYTHON_VER}"
echo "  EDTA:          ${EDTA_VER}"
echo "  验证报告:      ${VALIDATION}"
echo "  源码备份:      ${BACKUP_DIR}/"
echo ""
echo "  下一步:"
echo "    1. cat ${VALIDATION}     # 检查所有 ✓ ✗"
echo "    2. 核对必需工具和依赖是否全部可用"
echo "    3. sbatch 01_prepare_task_table.sh   # 准备全部样本/染色体任务表"
echo ""
echo "  激活方式: conda activate ${ENV_NAME}"
echo "============================================================"
