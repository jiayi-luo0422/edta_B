#!/bin/bash
#SBATCH --job-name=panEDTA_final
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --qos=normal
#SBATCH --output=logs/panEDTA_final_%j.out
#SBATCH --error=logs/panEDTA_final_%j.err

#==========================================================================
# 05_build_pan_te_library.sh   (正式版: 固定样本表)  ——【官方复用模式】
#
# 关键设计 (取代软链伪装版, 杜绝重跑 de novo):
#   本版 panEDTA(v2.2.2) 对 genome 做 realpath(会解软链); 它用
#     <G>.mod.EDTA.TEanno.sum           作为"已有注释→跳过 de novo"的触发文件
#     <G>.mod.EDTA.anno/<G>.mod.EDTA.RM.out  喂 find_flTE 建 pan 库
#   因此本脚本:
#     1) 基因组以【硬链】摆进 stage/ (规避 realpath 把路径带出 stage);
#     2) 按 EDTA 原生命名预置三件套:
#          <G>.mod.EDTA.TEanno.sum                 非空占位, 仅触发, 内容不解析
#          <G>.mod.EDTA.TElib.fa                   → merged 单品种 reclassified 库
#          <G>.mod.EDTA.anno/<G>.mod.EDTA.RM.out   → repeatmasker_self .out 结果
#     3) panEDTA 以 -a 0 只建 pan 库 (跳过 Step3 整基因组重注释——那步走
#        --step final 需 EDTA raw/combine 中间件, per-chr 工作流没有;
#        重注释等 pan 库定稿后由06统一跑 RepeatMasker)。
#   若复用判据未命中(日志出现 De novo annotate)或未产出库, 自动回退 cd-hit+mmseqs。
#==========================================================================

set -euo pipefail
SCRIPT_DIR=$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)
source "${SCRIPT_DIR}/workflow_common.sh"
activate_edta

THREADS=32
WORK_DIR="${PANTE_DIR}/panEDTA_work"                # ★ 正式工作目录

LIB_SUFFIX="_TElib_nr.reclassified.fa"

#--------------------------------------------------------------------------
# 从固定样本表读取，与02–04使用相同批次。
#--------------------------------------------------------------------------
SAMPLES=()
while IFS= read -r sample || [ -n "$sample" ]; do
    [ -z "$sample" ] || SAMPLES+=("$sample")
done < "${SAMPLE_TABLE}"

if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "✗ 样本表为空: ${SAMPLE_TABLE}, 退出"; exit 1
fi

mkdir -p "${WORK_DIR}"
# 本次工具输出独立存放，失败不能误用上一轮残留文件。
WORK_DIR=$(mktemp -d "${WORK_DIR}/run.XXXXXX")
cd "${WORK_DIR}"

echo "============================================================"
echo "[$(date)] Pan-TE: panEDTA 复用模式 (正式版, ${#SAMPLES[@]} 个品种)"
echo "  品种:        ${SAMPLES[*]}"
echo "  WORK_DIR:    ${WORK_DIR}"
echo "  THREADS:     ${THREADS}"
echo "  lib suffix:  ${LIB_SUFFIX}"
echo "============================================================"

########################################################################
# 1. 找 panEDTA.sh
########################################################################
if EDTA_BIN=$(command -v EDTA.pl); then
    EDTA_PREFIX=$(dirname "${EDTA_BIN}")/..
    EDTA_PREFIX=$(realpath "${EDTA_PREFIX}")
    EDTA_BASE="${EDTA_PREFIX}/share/EDTA"
    PANEDTA_SH="${EDTA_BASE}/panEDTA.sh"
    
    if [ ! -f "${PANEDTA_SH}" ]; then
        echo "  ⚠ panEDTA.sh 不在 ${PANEDTA_SH}, 搜索其它位置..."
        PANEDTA_SH=$(find "${EDTA_PREFIX}" -name "panEDTA.sh" -type f 2>/dev/null | head -1)
        if [ -z "${PANEDTA_SH}" ]; then
            echo "  ✗ panEDTA.sh 未找到, 回退自定义合并"; USE_CUSTOM=1
        else
            echo "  ✓ 找到: ${PANEDTA_SH}"; USE_CUSTOM=0
        fi
    else
        echo "  ✓ panEDTA.sh: ${PANEDTA_SH}"; USE_CUSTOM=0
    fi

else
    echo "  ✗ EDTA.pl 未找到, 回退自定义合并"
    USE_CUSTOM=1
fi

########################################################################
# 2. 前置检查 + 准备输入清单
########################################################################
echo ""
echo "[$(date)] 前置检查"

n_total=${#SAMPLES[@]}
n_genome_ok=0; n_lib_ok=0; n_rm_ok=0

GENOME_LIST="${WORK_DIR}/genome_list.txt"
CDS_LIST="${WORK_DIR}/cds_list.txt"
TELIB_LIST="${WORK_DIR}/telib_list.txt"
RMOUT_LIST="${WORK_DIR}/rmout_list.txt"
> "${GENOME_LIST}"; > "${CDS_LIST}"; > "${TELIB_LIST}"; > "${RMOUT_LIST}"

for sample in "${SAMPLES[@]}"; do
    GENOME="${GENOME_DIR}/${sample}.chrom.fasta"
    TELIB="${EDTA_OUT_DIR}/${sample}/merged/${sample}${LIB_SUFFIX}"
    # RM.out 已移动至 repeatmasker_self 并重命名为 .out
    RMOUT="${EDTA_OUT_DIR}/${sample}/repeatmasker_self/${sample}.chrom.fasta.out"
    CDS="${CDS_DIR}/${sample}.cds.fa"
    RM_MARK="${EDTA_OUT_DIR}/${sample}/repeatmasker_self/repeatmasker.complete"

    status="${sample}:"
    if [ -s "${GENOME}" ]; then n_genome_ok=$((n_genome_ok+1)); status="${status} fa✓"; else status="${status} fa✗"; fi
    if current_reclassified_library "${sample}"; then
        n_lib_ok=$((n_lib_ok+1)); status="${status} lib✓"
    else
        status="${status} lib✗(未验收或上游已更新)"
    fi
    if [ -s "${RMOUT}" ] && [ -f "${RM_MARK}" ] &&
       [ ! "${GENOME}" -nt "${RM_MARK}" ] && [ ! "${TELIB}" -nt "${RM_MARK}" ] &&
       [ ! "${GENOME}" -nt "${RMOUT}" ] && [ ! "${TELIB}" -nt "${RMOUT}" ]; then
        n_rm_ok=$((n_rm_ok+1)); status="${status} rm✓"
    else
        status="${status} rm✗(未验收或比输入旧)"
    fi
    echo "  ${status}"

    echo "${GENOME}"  >> "${GENOME_LIST}"
    if [ -s "${CDS}" ]; then echo "${CDS}" >> "${CDS_LIST}"; fi
    if [ -s "${TELIB}" ]; then echo "${TELIB}" >> "${TELIB_LIST}"; fi
    if [ -s "${RMOUT}" ]; then echo "${RMOUT}" >> "${RMOUT_LIST}"; fi
done

echo ""
echo "  汇总: ${n_genome_ok}/${n_total} 基因组, ${n_lib_ok}/${n_total} 库, ${n_rm_ok}/${n_total} RM输出"
if [ "${n_genome_ok}" -lt "${n_total}" ] || [ "${n_lib_ok}" -lt "${n_total}" ] || [ "${n_rm_ok}" -lt "${n_total}" ]; then
    echo "  ✗ 数据不全, 退出 (检查样本名 / 03 重分类 / repeatmasker_self 目录及文件是否就绪)"
    echo "    上面带 ✗ 的品种就是缺料的, 补齐再投"
    exit 1
fi

########################################################################
# 3. 官方 panEDTA: 复用单品种注释, 仅建 pan 库 (-a 0), 绝不重跑 de novo
########################################################################
if [ "${USE_CUSTOM}" -eq 0 ]; then
    echo ""
    echo "[$(date)] 官方 panEDTA (复用模式, -a 0)"

    if [ ! -s "${CDS_LIST}" ]; then
        # 沿用本流程的分支条件：仅存在非空 CDS 时尝试 panEDTA。
        echo "  ⚠ 本流程无可用 CDS，沿用自定义合并分支"
        USE_CUSTOM=1
    else
        REF_CDS=$(head -1 "${CDS_LIST}")
        STAGE="${WORK_DIR}/stage"
        PANRUN="${WORK_DIR}/panrun"          # 实际跑目录, 必须 ≠ STAGE
        rm -rf "${STAGE}" "${PANRUN}"; mkdir -p "${STAGE}" "${PANRUN}"

        for sample in "${SAMPLES[@]}"; do
            G_REAL="${GENOME_DIR}/${sample}.chrom.fasta"
            TELIB="${EDTA_OUT_DIR}/${sample}/merged/${sample}${LIB_SUFFIX}"
            RMOUT="${EDTA_OUT_DIR}/${sample}/repeatmasker_self/${sample}.chrom.fasta.out"
            G="${sample}.chrom.fasta"        # = realpath 后 panEDTA 认的 <genome>

            # 基因组用硬链(规避 realpath 把路径带出 STAGE); 跨 FS 退化为 cp
            ln "${G_REAL}" "${STAGE}/${G}" 2>/dev/null || cp "${G_REAL}" "${STAGE}/${G}"
            # 触发"复用已有注释"判据(非空占位, panEDTA 不解析其内容)
            echo "reuse per-cultivar EDTA annotation" > "${STAGE}/${G}.mod.EDTA.TEanno.sum"
            # 库 + RM.out 软链成 EDTA 原生命名
            ln -sf "${TELIB}" "${STAGE}/${G}.mod.EDTA.TElib.fa"
            mkdir -p "${STAGE}/${G}.mod.EDTA.anno"
            ln -sf "${RMOUT}" "${STAGE}/${G}.mod.EDTA.anno/${G}.mod.EDTA.RM.out"
        done

        GLIST="${PANRUN}/genomes.list"; > "${GLIST}"
        for sample in "${SAMPLES[@]}"; do
            echo "${STAGE}/${sample}.chrom.fasta" >> "${GLIST}"
        done

        cd "${PANRUN}"
        echo "  cmd: bash ${PANEDTA_SH} -g ${GLIST} -t ${THREADS} -c ${REF_CDS} -a 0 -o 0"
        # v2.2.2 的 Step 1 允许可选软链失败；仅对 Step 2 建库启用失败传播。
        CHECKED_PANEDTA="${WORK_DIR}/panEDTA.checked.sh"
        if awk '
            /^path=\$\(dirname "\$0"\)/ {print "path=$PANEDTA_BASE"; path_found++; next}
            /^## Step 2,/ {print "set -e -o pipefail"; step_found++}
            {print}
            END {if (path_found != 1 || step_found != 1) exit 1}
        ' "${PANEDTA_SH}" > "${CHECKED_PANEDTA}"; then
            set +e
            PANEDTA_BASE=$(dirname "${PANEDTA_SH}") bash "${CHECKED_PANEDTA}" \
                -g "${GLIST}" -t "${THREADS}" -c "${REF_CDS}" -a 0 -o 0 \
                2>&1 | tee "${WORK_DIR}/panEDTA.log"
            PAN_STATUS=("${PIPESTATUS[@]}")
            set -e
            PAN_EXIT=${PAN_STATUS[0]}
            [ "${PAN_STATUS[1]}" -eq 0 ] || PAN_EXIT=${PAN_STATUS[1]}
        else
            echo "  ✗ 无法定位 panEDTA 建库步骤，回退自定义合并" | tee "${WORK_DIR}/panEDTA.log"
            PAN_EXIT=1
        fi
        echo ""
        echo "[$(date)] panEDTA 退出码: ${PAN_EXIT}"

        n_reuse=$(count_matches "Existing EDTA annotation found" "${WORK_DIR}/panEDTA.log")
        if [ "${PAN_EXIT}" -ne 0 ] || [ "${n_reuse}" -ne "${n_total}" ]; then
            echo "  ✗ 官方退出码或复用数量不合格 (${PAN_EXIT}, ${n_reuse}/${n_total}), 回退自定义合并"
            USE_CUSTOM=1
        fi
        # ★ 硬校验: 绝不应回到 de novo
        if grep -q "De novo annotate" "${WORK_DIR}/panEDTA.log"; then
            echo "  ✗ 仍触发 de novo! 复用判据未命中, 回退自定义合并"
            USE_CUSTOM=1
        fi

        PAN_LIB=$(find "${PANRUN}" -maxdepth 1 -name "*.panEDTA.TElib.fa" -size +1k 2>/dev/null | head -1)
        if [ -n "${PAN_LIB}" ] && [ -s "${PAN_LIB:-/nonexistent}" ] && [ "${USE_CUSTOM}" -eq 0 ]; then
            echo "  ✓ panEDTA 库: ${PAN_LIB}"
            echo "  复用命中: ${n_reuse}/${n_total} (应等于 ${n_total})"
        else
            echo "  ✗ 未产出 panEDTA.TElib.fa, 回退自定义合并"
            USE_CUSTOM=1
        fi
        cd "${WORK_DIR}"
    fi
fi

########################################################################
# 4. 回退方案: 自定义跨品种合并 (cd-hit + mmseqs 两步)
########################################################################
if [ "${USE_CUSTOM}" -eq 1 ]; then
    echo ""
    echo "[$(date)] 回退: 自定义合并 (cd-hit + mmseqs 两步)"
    cd "${WORK_DIR}"

    COMBINED="${WORK_DIR}/all_samples_TElib.fa"; > "${COMBINED}"
    for sample in "${SAMPLES[@]}"; do
        TELIB="${EDTA_OUT_DIR}/${sample}/merged/${sample}${LIB_SUFFIX}"
        sed "s/^>/>${sample}__/" "${TELIB}" >> "${COMBINED}"
    done
    N_RAW=$(count_matches "^>" "${COMBINED}")
    [ "${N_RAW}" -gt 0 ] || { echo "  ✗ 合并库无FASTA记录"; exit 1; }
    echo "  跨品种合并: ${N_RAW} 条"

    CDHIT_OUT="${WORK_DIR}/all_samples_cdhit.fa"
    echo "  [$(date)] cd-hit-est (-G 1 -aS 0.8 -aL 0.8)"
    cd-hit-est -i "${COMBINED}" -o "${CDHIT_OUT}" \
        -c 0.8 -G 1 -aS 0.8 -aL 0.8 -n 5 -d 0 -M 65000 -T ${THREADS} \
        > "${WORK_DIR}/cdhit_pan.log" 2>&1 || { echo "  ✗ cd-hit失败: ${WORK_DIR}/cdhit_pan.log"; exit 1; }
    N_CDHIT=$(count_matches "^>" "${CDHIT_OUT}")
    [ "${N_CDHIT}" -gt 0 ] || { echo "  ✗ cd-hit输出无FASTA记录"; exit 1; }
    echo "    ${N_RAW} → ${N_CDHIT}"

    LINCLUST="${WORK_DIR}/mmseqs_linclust"
    LINCLUST_REP="${LINCLUST}_rep_seq.fasta"
    rm -rf "${WORK_DIR}/linclust_tmp"; mkdir "${WORK_DIR}/linclust_tmp"
    echo "  [$(date)] mmseqs linclust (90%)"
    mmseqs easy-linclust "${CDHIT_OUT}" "${LINCLUST}" "${WORK_DIR}/linclust_tmp" \
        --min-seq-id 0.90 -c 0.80 --cov-mode 0 --threads ${THREADS} \
        > "${WORK_DIR}/linclust.log" 2>&1 || { echo "  ✗ mmseqs linclust失败: ${WORK_DIR}/linclust.log"; exit 1; }
    rm -f "${LINCLUST}_all_seqs.fasta"; rm -rf "${WORK_DIR}/linclust_tmp"
    N_LIN=$(count_matches "^>" "${LINCLUST_REP}")
    [ "${N_LIN}" -gt 0 ] || { echo "  ✗ linclust输出无FASTA记录"; exit 1; }
    echo "    ${N_CDHIT} → ${N_LIN}"

    CLUSTER="${WORK_DIR}/mmseqs_cluster"
    CLUSTER_REP="${CLUSTER}_rep_seq.fasta"
    rm -rf "${WORK_DIR}/cluster_tmp"; mkdir "${WORK_DIR}/cluster_tmp"
    echo "  [$(date)] mmseqs cluster (80%, sensitive)"
    mmseqs easy-cluster "${LINCLUST_REP}" "${CLUSTER}" "${WORK_DIR}/cluster_tmp" \
        --min-seq-id 0.80 -c 0.80 --cov-mode 0 \
        --cluster-mode 1 --cluster-reassign 1 --kmer-per-seq 80 -s 7.5 \
        --threads ${THREADS} > "${WORK_DIR}/cluster.log" 2>&1 || { echo "  ✗ mmseqs cluster失败: ${WORK_DIR}/cluster.log"; exit 1; }
    rm -f "${CLUSTER}_all_seqs.fasta"; rm -rf "${WORK_DIR}/cluster_tmp"

    PAN_LIB="${CLUSTER_REP}"
    N_FINAL=$(count_matches "^>" "${PAN_LIB}")
    [ "${N_FINAL}" -gt 0 ] || { echo "  ✗ cluster输出无FASTA记录"; exit 1; }
    echo "    ${N_LIN} → ${N_FINAL}"
    echo "  ⚠ 提示: 回退路把 family 压到 80%, 较激进; 若 TIP 发现同源插入被并到"
    echo "    一个 family, 考虑把这步阈值提到 0.90"
fi

########################################################################
# 5. 输出最终库 + 分类分布
########################################################################
if [ -z "${PAN_LIB:-}" ] || [ ! -s "${PAN_LIB}" ]; then
    echo "  ✗ 没有可用 pan 库 (官方与回退均失败), 退出"; exit 1
fi

FINAL_LIB="${PANTE_DIR}/Tea_panTE.lib.fasta"            # ★ 正式 pan 库
N_FINAL=$(count_matches "^>" "${PAN_LIB}")
[ "${N_FINAL}" -gt 0 ] || { echo "  ✗ pan库无FASTA记录"; exit 1; }
FINAL_NEW=$(mktemp "${PANTE_DIR}/.Tea_panTE.lib.XXXXXX")
cp "${PAN_LIB}" "${FINAL_NEW}"
mv "${FINAL_NEW}" "${FINAL_LIB}"

echo ""
echo "============================================================"
echo "[$(date)] Pan-TE 库完成 (正式版, ${n_total} 品种)"
echo "============================================================"
echo "  最终库:      ${FINAL_LIB}"
echo "  来源:        ${PAN_LIB}"
echo "  family 数:   ${N_FINAL}"

echo ""
echo "  分类分布:"
grep "^>" "${FINAL_LIB}" | awk -F'#' '{print $2}' | awk -F'/' '
{
    if($1=="LTR" && $2=="Gypsy")         gypsy++
    else if($1=="LTR" && $2=="Copia")    copia++
    else if($1=="LTR" && tolower($2)=="unknown")  ltrunk++
    else if($1=="LTR")                   ltr++
    else if($1=="LINE")                  line++
    else if($1=="SINE")                  sine++
    else if($1=="MITE")                  mite++
    else if($1=="DNA" && $2=="Helitron") helitron++
    else if($1=="DNA")                   dnatir++
    else if($1=="RC")                    rc++
    else                                 other++
}
END {
    printf "    LTR/Gypsy:     %5d\n", gypsy+0
    printf "    LTR/Copia:     %5d\n", copia+0
    printf "    LTR/unknown:   %5d\n", ltrunk+0
    printf "    LTR/其它:       %5d\n", ltr+0
    printf "    LINE:          %5d\n", line+0
    printf "    SINE:          %5d\n", sine+0
    printf "    MITE:          %5d\n", mite+0
    printf "    DNA/真TIR:      %5d\n", dnatir+0
    printf "    DNA/Helitron:  %5d\n", helitron+0
    printf "    RC:            %5d\n", rc+0
    printf "    其它/unknown:   %5d\n", other+0
}'

echo ""
echo "  ★ 完成判据:"
echo "    - 退出码 0 且 ${FINAL_LIB} 非空"
if [ "${USE_CUSTOM}" -eq 0 ]; then
    echo "    - 官方模式: panEDTA退出0, 复用 ${n_total}/${n_total}, 无 De novo annotate"
    echo "    - 各品种 flTE 贡献: wc -l ${WORK_DIR}/panrun/*.keep.list"
else
    echo "    - 回退模式: cd-hit与两步mmseqs均退出0，输出均已验收"
fi
echo "    - 峰值内存/耗时: sacct -j ${SLURM_JOB_ID:-JOB_ID} --format=MaxRSS,Elapsed"
echo ""
echo "  下一步: 见README第06步，提交06_mask_genomes_with_pan_te.sh统一重注释"
echo ""