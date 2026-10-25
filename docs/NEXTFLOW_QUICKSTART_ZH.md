# Katana 运行指南：先按手头的数据选入口

**可以从 bulk FASTQ、基因组 BAM、匹配的已有分析结果开始；已有结果直接复用。
完成下游后生成定稿七章报告，再按需选择 Gene/Event 做独立 deep dive。**

目前是复用 nf-core 模块的独立 Nextflow pipeline，尚不是获批的
`nf-core/katana` 正式发行版。这里按已发表 nf-core pipeline 的文档结构整理，
但参数与格式以 Katana 实际实现为准。

| 你已有的东西 | 用哪个入口 | 会做什么 |
| --- | --- | --- |
| bulk FASTQ | Nextflow `--mode complete` | 上游 → Katana 下游 → 七章报告。 |
| 基因组比对 BAM | 同一命令，换 BAM samplesheet | 保留已有比对，审核 BAM/参考，做所需定量、差异和下游报告。 |
| 只想生成上游产物 | `--mode upstream-only` | 生成完整原生输入 bundle 后停止。 |
| 已有完整 Katana 上游 bundle/project | 导入后 `--mode downstream-only` | 不重做读段处理，运行选定下游。 |
| 已完成下游结果 | `katana-flow report` / `katana-report` | 组装报告；不会补算缺失分析。 |
| 想深入某些基因/event | 浏览器勾选或 `katana-targets` | Gene card、Event card、Sashimi、独立报告、可导回浏览器的数据。 |

外部 `nf-core/rnaseq`、`rnasplice`、`differentialabundance` 的产物可以复用，
**需要准确绑定实际文件、样本、对比、参考和原方法**；不是把任意结果文件夹
交给报告命令就能补齐所有分析。rnaseq 本身也不产生完整 rMATS、网络和 MotifTrace。

## 1. 安装与最常用命令

先按[Release 安装说明](RELEASE_INSTALL.md)安装匹配的 wheel、加载容器，
将 `KATANA_IMAGE`、`DEXSEQ_IMAGE` 设为已加载/构建的真实镜像或 SIF。
Linux 执行端需要 Nextflow ≥25.10.4、兼容 Java 和对应引擎。软件依赖由容器提供；
FASTA/GTF、motif、GO/KEGG、CLIP 等资源仍由项目明确配置。

```sh
katana-flow export-pipeline --outdir tools/katana-pipeline

nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode complete --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --reference_build YOUR_BUILD --annotation_release YOUR_RELEASE \
  --project_template inputs/project-template.json --flow_config inputs/flow.json \
  --katana_container "$KATANA_IMAGE" --dexseq_tools_container "$DEXSEQ_IMAGE" \
  --outdir results -work-dir work
```

FASTQ/BAM 用同一条入口命令，由表中字段决定。配置、输入目录与 `results/work`
分开放；以上命令在研究项目根目录运行，不包含任何本机或 Dardel 固定路径。
完整说明见 [Nextflow usage](NEXTFLOW_USAGE.md)。

**主分支源码中的自动判断链特异性：** `strandedness=auto` 请使用当前
`main` 源码中的 `pipeline/main.nf`：

```sh
git clone --branch main https://github.com/Taratype0/Katana-splicing-tool.git katana-source
nextflow run katana-source/pipeline/main.nf --help
```

随后将上面运行命令的 `tools/katana-pipeline/main.nf` 换成
`katana-source/pipeline/main.nf`，保留其他真实输入和匹配环境参数。
已发布的 `v2.0.0-beta.5.dev20261025` wheel 和原有导出包保持冻结，仍要求手动填写
链特异性；主分支源码更新不会改变这些已发布文件。原有下游容器仍可使用；
新增判断步骤使用独立固定版本的 RSeQC/gffread 镜像。

## 2. 样本表怎么填

可直接复制：[FASTQ 表](../examples/nextflow/samples_fastq.csv)、
[BAM 表](../examples/nextflow/samples_bam.csv)、
[简单对比表](../examples/nextflow/contrasts.csv)、
[12 个指定对比示例](../examples/nextflow/contrasts_factorial.csv)。
请先阅读[模板说明](../examples/nextflow/README.md)，替换真实文件路径和样本身份。

FASTQ 示例，实际每组写齐全部独立生物学样本：

```csv
sample,condition,fastq_1,fastq_2,strandedness,read_length,subject
Control_1,Control,../reads/Control_1_R1.fastq.gz,../reads/Control_1_R2.fastq.gz,reverse,150,animal_01
Treated_1,Treated,../reads/Treated_1_R1.fastq.gz,../reads/Treated_1_R2.fastq.gz,reverse,150,animal_04
```

BAM 示例：

```csv
sample,condition,genome_bam,strandedness,single_end,read_length,subject
Control_1,Control,../alignments/Control_1.bam,reverse,false,150,animal_01
Treated_1,Treated,../alignments/Treated_1.bam,reverse,false,150,animal_04
```

- `sample` 唯一且以字母开头，只用字母、数字、下划线、点、连字符。
- `condition` 以字母开头，只用字母、数字、下划线；与 contrast 表完全一致。
- FASTQ 必须 gzip 压缩；`fastq_2` 留空表示单端，主验收覆盖双端。
- `strandedness` 填 `forward/reverse/unstranded`；上述 `main` 源码入口也接受 `auto`。
  手动指定值保持不变，可与 `auto` 行混用。
- BAM 必须是基因组坐标的 BAM，明确 `single_end`；转录组 BAM 不能替代。
- `read_length` 为名义读长，可用 `--read_length` 补充；目前整套研究须一致。
- `subject` 是真实生物学单位；独立动物填不同 ID，不能把 lane 当生物学重复。
- 读段路径相对于 **samplesheet 所在目录**，因此 `inputs/samples.csv` 用 `../reads/`。

`auto` 对 FASTQ、基因组 BAM 使用同一条判断路线：FASTQ 先经 STAR 对齐，现成 BAM
不重新对齐；共享 BAM/参考审计、排序、索引后，由 gffread 将匹配的**完整未过滤 GTF**
转成 BED12，再由 RSeQC 5.0.3 `infer_experiment.py` 判断。仅 `auto` 样本增加这一步。
默认读取最多 200,000 条合格、与基因模型重叠的比对记录，MAPQ 至少 30；这是按 BAM
顺序选取的读段记录，并非随机抽样或 200,000 个双端片段。

正向、反向比例先在可判断记录中归一化：一侧至少 0.8 判为该方向；否则两侧比例差
不超过 0.1 判为 `unstranded`。其他情况、没有可判断证据、格式/单双端不符均报错，
**不会默认改成 unstranded**。所有所选样本判断完成、每个所选对比实际链特异性一致后，
才进入 Salmon、DEXSeq、rMATS。请求值、实际值、原始 RSeQC 文本/日志及失败诊断
保存在 pipeline QC 和 bundle provenance 中。此步骤不能还原 BAM 上游已丢弃的读段。

**当前原始数据入口要求每个样本一行、一套已合并文件。** 不支持像 rnaseq 那样
用重复 sample 行自动合并多 lane；先合并相容 lane，保留映射记录，或接入已合并 BAM。
元数据检查支持 lane 身份盘点，不等于原始数据 DAG 已实现 lane 合并。

上游差异检验至少每组 2 个独立样本；完整报告默认逐对比建网，需每个对比至少
6 个独立样本。因此 2+2 可以达到上游最低要求，却不能通过默认完整下游。

## 3. Contrast 表格式与选择

```csv
contrast,treatment,control
Treated_vs_Control,Treated,Control
```

只填你希望比较的行，不穷尽所有组合。`contrast` 是唯一名称，真正的方向由
`treatment/control` 决定，不能靠名字猜。若有四种背景与两个处理水平，可选：

| 对比类型 | 要解决的数值问题 |
| --- | --- |
| 各处理下 `KO vs CONTROL` | 相同处理下两种背景的差异。 |
| 各处理下 `WTrescue vs KO`、`MUTrescue vs KO` | 两种指定序列相对 KO 的差异。 |
| 各处理下 `MUTrescue vs WTrescue` | 两种指定序列之间的差异。 |
| 各背景下 `DP vs NO_DP` | 该背景内的处理响应。 |

这样是 8 个背景差异＋4 个处理响应，共 12 个选定比较；可按项目删减。
[完整可复制 CSV](NEXTFLOW_USAGE.md#4-contrastsheet-input-and-direction)。
软件不会从 `rescue` 名字自行判断是否恢复功能，也不会把响应之差伪装成交互模型。

原始数据完整路线目前只支持独立样本 `~ condition`，**不支持**自动配对、批次校正
或任意交互公式。`subject/batch` 可记录信息，但不会默默加入统计模型。
可选 `paired/blocking/design/formula` 列中的不支持模型会被拒绝。
外部已拟合结果保留原模型；专门的下游交互分析使用独立显式设计。

DESeq2 为 `log2(treatment/control)`；rMATS 为 treatment/group1 减 control/group2。
SUPPA 由适配器固定输入顺序。DEXSeq 以实际原始 fold-change 列名为准，当前原生导出
可能是反方向，不能统一强行改符号。JC/JCEC 与每组样本向量顺序必须明确保留。

## 4. 哪些参数需要选择

| 选择 | 当前默认 / 配置位置 |
| --- | --- |
| 做哪些对比 | 上游 `contrasts.csv`；下游 `flow.json:contrast_ids` 再选已声明子集。 |
| 参考版本 | `--fasta/--gtf/--reference_build/--annotation_release`。`--genome` 只是备用标签，不自动下载。 |
| 自动链特异性，源码入口 `auto` 行 | `--stranded_threshold=0.8`、`--unstranded_threshold=0.1`；RSeQC `--strandedness_sample_size=200000`、`--strandedness_mapq=30`。不改变手动指定行。 |
| 排除假基因、线粒体编码基因 | 原生分析参考中默认执行，先于定量/差异；完整原注释和排除记录保留。 |
| 普通 DAS / RS-DAS | `flow.json:event_set`，默认 `stable_significant` 即 RS-DAS；普通为 `significant`。 |
| 读数支持 | 固定规则：选定对比中**每个样本** `IJC+SJC≥10`；不以平均值代替，无旧版额外一致性/近似 SE 门槛。 |
| DEG 阈值 | `parameters.differential.deg`：`fdr=0.05`、`min_abs_log2fc=1.0`，可自定义。 |
| DAS 显著集合 | `parameters.quality.stable_as`：`fdr=0.05`、`min_abs_dpsi=0.10`，可调；`differential.as` 继承已完成的 Quality 结果，不是 top-level `--das_fdr` 参数。 |
| rMATS 原生检验 | `--rmats_splice_diff_cutoff=0.0001` 等；这个参数不是下游 ΔPSI 0.10 筛选。 |
| JC/JCEC | 模板 `rmats_mode`；网络可明确 `network.parameters.rmats_mode`。记录真实模式。 |
| 网络 | `network.scope=contrast`，每个对比独立建网；模块最小 20 个事件、合并高度 0.25，可自定义。 |
| Bootstrap / Hub | 默认不做 bootstrap；`hub_metric=abs_kME`。连续指标和完整排序保留，可选标签不硬筛。 |
| MotifTrace | 每对比独立、多个 motif；正文最多 10 个显著模型及序列、P/FDR，完整模型和命中基因在输出里。 |
| GO/KEGG 背景 | 与该证据层相匹配的检测/可检验基因，不是所有注释基因；输入、背景和命中基因有下载入口。 |

阈值依据见 [THRESHOLDS.md](THRESHOLDS.md)；
**完整参数、默认值、上下游分别在哪里改**见 [NEXTFLOW_PARAMETERS.md](NEXTFLOW_PARAMETERS.md)。
完整报告依赖真实资源；Conda 安装不能替代缺失的 CLIP、motif 数据库或参考。

先检查再提交：

```sh
python tools/katana-pipeline/bin/katana_validate_inputs.py \
  --samplesheet inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf --outdir checks/raw
katana-flow preflight --project-template inputs/project-template.json \
  --config inputs/flow.json --out checks/readiness.json
```

## 5. 已有结果怎么接

```sh
katana-flow import --bundle upstream-results/katana_bundle \
  --project-template inputs/project-template.json --outdir imported-project

nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode downstream-only --project imported-project/project.json \
  --flow_config inputs/flow.json --katana_container "$KATANA_IMAGE" \
  --outdir downstream-results -work-dir downstream-work
```

`import` 接受完整、带 SHA 清单的 Katana 原生 bundle，不直接接受任意 nf-core 目录。
外部结果先明确映射到兼容 project/source；历史结果的只读采用、身份核验和缺项补算
见 [EXISTING_RESULTS.md](EXISTING_RESULTS.md)。

只生成报告，不运行统计计算：

```sh
katana-flow report --collection completed-collection --outdir report-final
katana-report --project results.json --outdir report-from-results \
  --compact --archive downloads/study-report.zip
```

第一条要求已完成 collection；第二条的 `results.json` 是明确绑定已完成结果的描述文件。
输出保留定稿七章和对比选择器。HTML 出现不等于完整验收；严格完整报告应有成功的
`flow.acceptance.json`，失败时保留诊断而不冒充完整结果。

## 6. Gene/Event deep dive 与结果导回

浏览器勾选目标 → 导出清单 ZIP → 本地/服务器绑定可访问 BAM 和索引 → 运行：

```sh
katana-targets inspect --request targets.zip
katana-targets plan --request targets.zip --bindings bindings.json
katana-targets run --request targets.zip --bindings bindings.json --outdir target-results
katana-targets results --outdir target-results
```

生成独立 `INDEX.html`、逐目标报告、Gene/Event cards、已有网络/motif 证据和实际
Sashimi（若请求且资源齐全），以及完整 TSV/Parquet/JSON、ZIP。不重做差异分析。
BAM 留在执行机器，不下载进报告包。

取回结果后打开 HTML 或在 Katana 浏览器导入；需要继续 CLI 工作时：

```sh
katana-targets restore --results completed-targets.zip --outdir restored-targets
katana-targets results --outdir restored-targets
katana-open
```

`bindings.json` 用 inspect 返回的实际资源 ID 绑定本机路径。
[完整目标运行与 Python API](TARGET_COMMANDS.md)。离线 HTML 能看结果；运行新 Sashimi
需要 Python 服务/CLI 和可访问 BAM，不能仅靠静态网页执行。

## 7. 其他运行方法、并行和恢复

- **Docker**：统一 Nextflow 路线，软件环境在镜像内。
- **Slurm＋Apptainer/Singularity**：`-profile slurm,apptainer`，配置真实账号/分区和节点可见镜像，需站点适配。
- **独立 Conda**：只想使用下游时，安装锁定的 Linux 原生环境，运行相同 CLI/API。
- **浏览器**：`katana-open`；服务器用 `katana-open --no-open`，通过安全隧道连接。
- **单模块 CLI/Python API**：按已配置项目运行某个模块，不强制整条重算。
- **目标任务批处理**：`katana-targets tasks --request targets.zip --format tsv`，在任意调度器任务中按真实 target ID 执行，每任务独立输出目录。

逐对比网络、MotifTrace、MAPP 按 DAG 依赖独立调度；Slurm 决定实际并发，没有隐藏
的 4 任务上限。`--max_cpus/max_memory/max_time` 是单任务资源上限，不是并发数。
不要在共享登录节点用 local executor 跑重计算。

原命令加 `-resume` 复用匹配缓存；保留原 `.nextflow`、work、输入和已执行 plan 原路径。
换镜像/代码/资源/参数可能改变对应任务身份，不应改缓存强迫命中。清理 work 前先确认
已验证的独立最终结果齐全；清理后将失去对应断点续跑能力。

参考结构：[rnaseq](https://nf-co.re/rnaseq/usage)、
[rnasplice](https://nf-co.re/rnasplice/usage)、
[differentialabundance](https://nf-co.re/differentialabundance/usage)。
这些官方 pipeline 的全部参数不能原样套用到 Katana；这里记录的是当前实际可执行接口。
