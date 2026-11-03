# Edit project paths here before preparing task tables or submitting jobs.
# EDTA_CONFIG can point to another copy of this file.
PROJDIR=${PROJDIR:-/public/home/a326/luojy/project/graph_pangenome_LJ43_hap1}
GENOME_DIR=${GENOME_DIR:-$PROJDIR/data/chromosome_only}
EDTA_OUT_DIR=${EDTA_OUT_DIR:-$PROJDIR/data/edta_new}
CDS_DIR=${CDS_DIR:-$PROJDIR/cds}
PANTE_DIR=${PANTE_DIR:-$PROJDIR/04_panTE}
PANTE_MASK_DIR=${PANTE_MASK_DIR:-$PROJDIR/05_panTE_annotation}
TASK_TABLE=${TASK_TABLE:-$EDTA_OUT_DIR/edta_tasks.tsv}
SAMPLE_TABLE=${SAMPLE_TABLE:-$EDTA_OUT_DIR/edta_samples.txt}
CONDA_SH=${CONDA_SH:-/public/home/a326/luojy/miniconda3/etc/profile.d/conda.sh}
SAMTOOLS_BIN=${SAMTOOLS_BIN:-/public/home/a326/luojy/miniconda3/envs/samtools_env/bin/samtools}
EDTA_ENV=${EDTA_ENV:-EDTA2_fix}
TESORTER_ENV=${TESORTER_ENV:-tesorter}
