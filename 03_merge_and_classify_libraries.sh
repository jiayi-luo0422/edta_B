#!/bin/bash
#SBATCH --job-name=merge_reclass_v2
#SBATCH --partition=core128
#SBATCH --account=a326
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --array=1
#SBATCH --qos=normal
#SBATCH --output=logs/merge_reclass_%A_%a.out
#SBATCH --error=logs/merge_reclass_%A_%a.err

#==========================================================================
# 03_merge_and_classify_libraries.sh
#
# 单品种一条龙: 跨 chr 合并 -> cd-hit-est 去冗余 -> TEsorter 重分类 LTR/unknown
#   产物(最终给 04 用): 02_te_annotation_v2/<sample>/merged/<sample>_TElib_nr.reclassified.fa
#
# 两个独立断点, 各自幂等:
#   Stage A 守 nr.complete + <sample>_TElib_nr.fa            (合并 + 去冗余, 中间库)
#   Stage B 守 reclass.complete + <sample>_TElib_nr.reclassified.fa (重分类后, 最终库)
#   => 标记与输入时间戳都满足才复用；无标记旧库在本次成功前保留。
#
# 前置: 样本表中本品种全部chr的02结果已验收 (TElib.fa + edta.complete)
# 环境: EDTA2_fix(cd-hit) + tesorter(TEsorter), 按阶段切换
# 全量投递: 按 edta_samples.txt 行数 N 使用 --array=1-N%4
# 完成后: 按 edta_samples.txt 行数 N 提交 sbatch --array=1-N%4 04_mask_genomes_with_repeatmasker.sh
#
# 设计要点:
#   - cd-hit 参数: -G 1 -aS 0.8 -aL 0.8 -c 0.8 (EDTA 源码标准)
#   - TEsorter: rexdb-plant, 只收 Gypsy/Copia 回填 class, 其余仍 unknown
#   - 不做预过滤(LTR 全保留, panEDTA 阶段处理)
#==========================================================================

set -euo pipefail
SCRIPT_DIR=$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)
source "${SCRIPT_DIR}/workflow_common.sh"

THREADS=16
TESORTER_DB="rexdb-plant"

########################################################################
# 1. 从固定样本表选品种
########################################################################
select_sample
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

PERCHR_DIR="${EDTA_OUT_DIR}/${SAMPLE}/per_chr"
MERGE_DIR="${EDTA_OUT_DIR}/${SAMPLE}/merged"
WORK="${EDTA_OUT_DIR}/${SAMPLE}/tesorter_reclass"
GENOME="${GENOME_DIR}/${SAMPLE}.chrom.fasta"
GENOME_FAI="${GENOME}.fai"

NR_OUT="${MERGE_DIR}/${SAMPLE}_TElib_nr.fa"                       # Stage A 产物(中间库)
RECLASS_OUT="${MERGE_DIR}/${SAMPLE}_TElib_nr.reclassified.fa"     # Stage B 产物(最终库)

NR_MARK="${MERGE_DIR}/nr.complete"
RECLASS_MARK="${MERGE_DIR}/reclass.complete"
mkdir -p "${MERGE_DIR}"
[ -s "${GENOME_FAI}" ] || { echo "  ✗ genome .fai 缺失: ${GENOME_FAI}"; exit 1; }
A_CURRENT=1
A_REBUILT=0
[ -s "${NR_OUT}" ] && [ -f "${NR_MARK}" ] || A_CURRENT=0
if [ "${GENOME_FAI}" -nt "${NR_MARK}" ] || [ "${NR_OUT}" -nt "${NR_MARK}" ]; then A_CURRENT=0; fi
while read -r chr rest; do
    TELIB="${PERCHR_DIR}/${chr}/${chr}.fa.mod.EDTA.TElib.fa"
    CHR_MARK="${PERCHR_DIR}/${chr}/edta.complete"
    if ! current_chromosome_result "$SAMPLE" "$chr"; then
        echo "  ✗ ${chr} 缺少已验收的02结果；请先核验/完成02 (edta.complete)"; exit 1
    fi
    if [ "${TELIB}" -nt "${NR_MARK}" ] || [ "${CHR_MARK}" -nt "${NR_MARK}" ]; then A_CURRENT=0; fi
done < "${GENOME_FAI}"

echo "============================================================"
echo "[$(date)] 单品种合并+重分类: ${SAMPLE} (task ${TASK_ID})"
echo "============================================================"

########################################################################
# Stage A: 跨 chr 合并 + cd-hit-est 去冗余  -> NR_OUT
########################################################################
if [ "${A_CURRENT}" -eq 1 ]; then
    echo "  [A 跳过] 已有去冗余库: $(grep -c '^>' "${NR_OUT}") families"
else
    activate_edta

    # --- A1. 检查 chr 完成度 ---
    echo ""
    echo "[$(date)] [A] 检查 chr 完成度"
    n_chr=0; n_done=0; n_failed=0; MISSING=""
    while read -r chr rest; do
        n_chr=$((n_chr + 1))
        TELIB="${PERCHR_DIR}/${chr}/${chr}.fa.mod.EDTA.TElib.fa"
        LOG="${PERCHR_DIR}/${chr}/edta_${chr}.log"
        if [ ! -s "${TELIB}" ]; then
            MISSING="${MISSING} ${chr}"; continue
        fi
        if [ -s "${LOG}" ]; then
            n_daemon=$(count_matches "daemonic processes are not allowed" "${LOG}")
            if [ "${n_daemon}" -gt 0 ]; then
                n_failed=$((n_failed + 1))
                echo "  ⚠ ${chr}: TElib 存在但 log 含 daemonic, 不合并"
                continue
            fi
        fi
        n_done=$((n_done + 1))
    done < "${GENOME_FAI}"

    echo "  ${SAMPLE}: ${n_done}/${n_chr} 干净, ${n_failed} 失败"
    [ -n "${MISSING}" ] && echo "  缺失 chr:${MISSING}"
    if [ "${n_done}" -lt "${n_chr}" ]; then
        echo "  ✗ ${SAMPLE} 未全部 chr 完成, 退出 (重投 03 补齐后再来)"
        exit 1
    fi

    # --- A2. 跨 chr 合并 (header 加 chr 前缀防 ID 冲突) ---
    echo ""
    echo "[$(date)] [A] 跨 chr 合并"
    A_WORK=$(mktemp -d "${MERGE_DIR}/.nr_run.XXXXXX")
    COMBINED="${A_WORK}/all_chr_TElib_combined.fa"
    NR_NEW="${A_WORK}/$(basename "${NR_OUT}")"
    > "${COMBINED}"
    while read -r chr rest; do
        TELIB="${PERCHR_DIR}/${chr}/${chr}.fa.mod.EDTA.TElib.fa"
        [ -s "${TELIB}" ] && sed "s/^>/>${chr}__/" "${TELIB}" >> "${COMBINED}"
    done < "${GENOME_FAI}"
    N_RAW=$(count_matches "^>" "${COMBINED}")
    [ "${N_RAW}" -gt 0 ] || { echo "  ✗ 合并库无FASTA记录"; exit 1; }
    echo "  合并后: ${N_RAW} 条"

    # --- A3. cd-hit-est 去冗余 ---
    echo ""
    echo "[$(date)] [A] cd-hit-est (-G 1 -aS 0.8 -aL 0.8 -c 0.8)"
    cd-hit-est \
        -i "${COMBINED}" -o "${NR_NEW}" \
        -c 0.8 -G 1 -aS 0.8 -aL 0.8 -n 5 -d 0 \
        -M 35000 -T ${THREADS} \
        > "${MERGE_DIR}/cdhit.log" 2>&1 || {
            echo "  ✗ cd-hit 失败，旧库保留"; tail -20 "${MERGE_DIR}/cdhit.log"; exit 1
        }

    if [ ! -s "${NR_NEW}" ]; then
        echo "  ✗ cd-hit 失败"; tail -20 "${MERGE_DIR}/cdhit.log"; exit 1
    fi
    N_NR=$(count_matches "^>" "${NR_NEW}")
    [ "${N_NR}" -gt 0 ] || { echo "  ✗ cd-hit输出无FASTA记录"; exit 1; }
    ratio=$(awk -v r=${N_RAW} -v n=${N_NR} 'BEGIN{printf "%.2f", r/n}')
    echo "  raw ${N_RAW} -> nr ${N_NR}  (压缩 ${ratio}×)"
    mv "${NR_NEW}" "${NR_OUT}"
    [ ! -f "${NR_NEW}.clstr" ] || mv "${NR_NEW}.clstr" "${NR_OUT}.clstr"
    printf "complete\n" > "${NR_MARK}"
    A_REBUILT=1
    rm -rf "${A_WORK}"
fi

########################################################################
# Stage B: TEsorter 重分类 LTR/unknown  -> RECLASS_OUT
########################################################################
if [ "${A_REBUILT}" -eq 0 ] && [ -s "${RECLASS_OUT}" ] && [ -f "${RECLASS_MARK}" ] &&
   [ ! "${NR_OUT}" -nt "${RECLASS_MARK}" ] && [ ! "${NR_MARK}" -nt "${RECLASS_MARK}" ] &&
   [ ! "${RECLASS_OUT}" -nt "${RECLASS_MARK}" ]; then
    echo "  [B 跳过] 已有重分类库: $(grep -c '^>' "${RECLASS_OUT}") families"
else
    activate_tesorter

    # 本次重分类独立目录，避免残留文件被误认为本次输出
    mkdir -p "${WORK}"
    B_WORK=$(mktemp -d "${WORK}/run.XXXXXX")
    cd "${B_WORK}"
    RECLASS_NEW="${B_WORK}/$(basename "${RECLASS_OUT}")"

    # --- B1. 抽 LTR/unknown consensus (去 #class 给 TEsorter) ---
    UNK="unknown_consensus.fa"
    awk '
    /^>/{ keep=(tolower($0) ~ /#ltr\/unknown([[:space:]]|$)/); if(keep){h=$0; sub(/#.*/,"",h); print h} next }
    keep { print }
    ' "${NR_OUT}" > "${UNK}"
    N_UNK=$(count_matches '^>' "${UNK}")
    echo "  [B] LTR/unknown family: ${N_UNK}"

    if [ "${N_UNK}" -eq 0 ]; then
        echo "  [B] 无 unknown, 直接复制原库为 reclassified (下游路径统一)"
        cp "${NR_OUT}" "${RECLASS_NEW}"
    else
        # --- B2. TEsorter ---
        TEsorter "${UNK}" -db "${TESORTER_DB}" -p "${THREADS}" || { echo "  ✗ TEsorter失败，旧重分类库保留"; exit 1; }
        CLS="${UNK}.${TESORTER_DB}.cls.tsv"
        [ -s "${CLS}" ] || { echo "  ✗ TEsorter 未出 cls.tsv"; exit 1; }

        # --- B3. name -> 新 class 映射 (只收 Gypsy/Copia) ---
        awk -F'\t' 'NR>1{
            if($3=="Gypsy")      print $1"\tLTR/Gypsy"
            else if($3=="Copia") print $1"\tLTR/Copia"
        }' "${CLS}" > reclass_map.tsv
        N_MAP=$(wc -l < reclass_map.tsv)
        echo "  [B] 解开(family): ${N_MAP} / ${N_UNK}  (仍 unknown: $((N_UNK - N_MAP)))"

        # --- B4. 改库 header -> RECLASS_OUT ---
        awk '
        BEGIN{ while((getline l < "reclass_map.tsv")>0){ split(l,a,"\t"); m[a[1]]=a[2] } }
        /^>/{ h=substr($0,2); split(h,p,"#"); if(p[1] in m) print ">"p[1]"#"m[p[1]]; else print $0; next }
        { print }
        ' "${NR_OUT}" > "${RECLASS_NEW}"
    fi

    [ -s "${RECLASS_NEW}" ] || { echo "  ✗ 重分类库写出失败"; exit 1; }
    [ "$(count_matches '^>' "${RECLASS_NEW}")" -eq "$(count_matches '^>' "${NR_OUT}")" ] || { echo "  ✗ 重分类family数不一致"; exit 1; }
    mv "${RECLASS_NEW}" "${RECLASS_OUT}"
    printf "complete\n" > "${RECLASS_MARK}"
    echo "  [B] ✓ 写出最终库: ${RECLASS_OUT}  ($(grep -c '^>' "${RECLASS_OUT}") families)"
fi

########################################################################
# 报告: 最终库分类分布
########################################################################
echo ""
echo "============================================================"
echo "  ${SAMPLE} 完成 (最终库: $(basename "${RECLASS_OUT}"))"
echo "------------------------------------------------------------"
grep "^>" "${RECLASS_OUT}" | awk -F'#' '{print $2}' | awk -F'/' '
{
    if($1=="LTR" && $2=="Gypsy")        gypsy++
    else if($1=="LTR" && $2=="Copia")   copia++
    else if($1=="LTR" && tolower($2)=="unknown") ltrunk++
    else if($1=="LTR")                  ltr++
    else if($1=="LINE")                 line++
    else if($1=="SINE")                 sine++
    else if($1=="MITE")                 mite++
    else if($1=="DNA" && $2=="Helitron") helitron++
    else if($1=="DNA")                  dnatir++
    else if($1=="RC")                   rc++
    else                                other++
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
echo "============================================================"
echo "[$(date)] ${SAMPLE} 全部完成 -> 可投 04_mask_genomes_with_repeatmasker.sh"
