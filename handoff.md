# 服务器B edta 01–06交接

有效代码基线：本地新版edta 01–06主线。通用版保留原集群配置，B适配版使用下列B路径；服务器实际脚本目录名称scripts4TE不代表旧版代码。

PROJDIR=/public/home/a326/luojy/project/graph_pangenome_LJ43_hap1；提交目录=$PROJDIR/scripts4TE；EDTA_CONFIG=$PROJDIR/scripts4TE/edta_config.sh。B默认core128、account a326、qos normal，日志位于scripts4TE/logs/。输入data/chromosome_only，任务表data/edta_new/edta_tasks.tsv，样本表同目录edta_samples.txt；02–04结果data/edta_new，05最终库04_panTE，06统一注释05_panTE_annotation。

已验收：10份FASTA及索引已传输，01为COMPLETED/0:0，10份样本、150项任务。Perl5.32.1/BioPerl1.7.8可用；SAMTOOLS_BIN=/public/home/a326/luojy/miniconda3/envs/samtools_env/bin/samtools，已在EDTA环境完成提取验证。Conda切换期间临时关闭nounset，恢复原状态并传播失败码。

EDTA_ENV实际位置/public/home/a326/luojy/.local/share/mamba/envs/EDTA2_fix，EDTA2.2.2、TIR-Learner4.08。B端源码与包记录确认TIR4使用grf_rs，grfmite-rs0.3.0与tirvish-rs0.1.1已安装；不补装旧grf-main、dask、swifter。EDTA自带依赖检查退出0，但不等于全流程兼容。EDTA内已有TEsorter1.5.1；03前仍需确认rexdb-plant可用性及TESORTER_ENV配置。

故障处理：21954的LTR_FINDER因入口旁缺少bin/cut.pl，产生只有注释的finder.combine.scn；21957复用了该无效文件，已取消。B已将share/EDTA/bin/LTR_FINDER_parallel/bin链接到share/LTR_FINDER_parallel/bin，cut.pl与convert_ltr_finder2.pl可读非空。按用户要求删除了仅AJBC_GCAv2/Chr01的未验收工作目录，未备份，再提交21958_1。最新日志显示13:55重新提取Chr01（238MB）、加载cds/AJBC_GCAv2.cds.fa并启动LTR候选检测，尚无整条染色体完成证据。

当前目标总并发6：21958_1运行时剩余2–150项并发5，退出后剩余数组提高到6；尚未收到剩余数组JobID，不视为已提交。不要重复提交第1项，不重跑已验收01。按实际sacct状态、退出码、日志和产物验收02后再推进03/04。当前代码仅标为运行验证中。

后续26份尚待迁移；05/06当前仅消费新10份样本表。36份整合需在新版框架另行准备有效输入，不能把26份追加到当前任务表而重跑。保留inputs人工批次清单及来源表。
