# MRQI AD Glioma evidence atlas

**RP2.4.1 · EDA91 · 2026-10-24 · exploratory research**

本仓库提供阿尔茨海默病与胶质瘤／历史GBM反向表达研究的完整基因级汇总、限定解释、图表及来源诊断。本版同步EDA90全文、Word和5处图–段映射修正，明确实际投稿集合及公开／本地范围；没有新增DE、P／q、MR拟合或生物学验证。

目前证据支持**来源、脑区及历史肿瘤–宿主区室限定的探索性反向表达，伴随蛋白冲突、缺失和批次敏感性**。它不确立人群保护、现代WHO2021 IDH-wildtype GBM对健康人群的独立确认、同细胞拮抗或因果治疗靶点。

## EDA91实际投稿集合

- [公开Word](docs/AD_Glioma_RP241_EDA90_PublicCommit_Submission_Aligned_Manuscript_2026-10-24.docx)与[全文Markdown](docs/CURRENT_MANUSCRIPT_WORKING_DRAFT.md)为同一104段内容；Word原字节及16页WPS验收绑定于EDA90，不再重拟合或重排。
- [实际55图清单](docs/EDA91_ACTUAL_FIGURE_SET.tsv)、[55图图例](docs/EDA91_CURRENT_SUBMISSION_LEGENDS.md)、[实际211表职责清单](docs/EDA91_ACTUAL_TABLE_SET.tsv)与[交付边界](docs/EDA91_ATTACHMENT_SCOPE.md)。
- 55图=12主图组件+43补图，41幅有公开像素、14幅仅本地。另2幅S1/S5历史背景保留，故仓库仍含43幅汇总图。原编号不变，不以P值或方向筛选。
- 本地211文件实装不等于211公开文件；此仓库仍只有35原汇总+2完整基因表派生映射。期刊尚未指定或提交。

## 阅读入口

- [完整81基因表](results/Complete81_GroupLevel_AllEvidence.tsv)与[原43基因表](results/Original43_GroupLevel_AllEvidence.tsv)：各1207列，所有基因、原方向及其余数值保留。
- [全文工作稿与29条参考文献](docs/CURRENT_MANUSCRIPT_WORKING_DRAFT.md)：Methods／Results／Discussion分开解释组织表达和已归档MR。
- [43幅原像素汇总图](docs/SELECTED_AGGREGATE_FIGURE_GALLERY.md)、[主图图例](docs/MAIN_FIGURE_LEGENDS.md)、[补图图例](docs/SUPPLEMENTARY_FIGURE_LEGENDS.md)。
- [全部211项实际披露清单](docs/All211_Actual_Payload_Index_And_Local_Disclosure.tsv)及[全部57幅图披露状态](docs/All57_Figure_Disclosure_Status.tsv)：索引项数不等于随包文件数。
- [17路线范围](reproduction/All17_Route_Scope.tsv)、[计算与发布范围](docs/PUBLICATION_AND_REPRODUCTION_SCOPE.md)、[来源及官方入口](docs/DATA_SOURCES.tsv)。

## 研究对象与主要结果

发现层整合GSE125583的218AD／70对照和GSE165595的15人／30库配对肿瘤–宿主资料。共同16068基因构成32136个方向性IUT假设，BH保留508、BY保留130。GSE159699的12AD／10Old形成81，GSE13276的5名成人历史GBM配对资料进一步形成43（28R1、15R2）。**508→81→43是顺序候选形成；81／81或43／43不是无偏复制率。**

| 证据层 | 当前可报告内容 | 必须保留的限制 |
| --- | --- | --- |
| 历史肿瘤–宿主 | GSE165595含11历史4级及4例3级；GSE13276仅5成人配对、12数组 | 宿主不是健康人群；未完成现代IDH-wildtype／WHO2021重分类 |
| AD脑区 | BA17主家族未支持；两个海马来源有有限CA1 CSRP1关联 | 平台尺度、RIN、脑区与供者差异不消失；不能平均不同尺度效应或q值 |
| 单核与肿瘤状态 | AD11／7供者六谱系；两tumour-only来源的状态定位 | 供者而非细胞是独立单位；状态定位不是疾病对照或同细胞反向机制 |
| 蛋白 | 既有AD／CPTAC与NBB公开汇总提供部分方向背景 | 44来源行、33／22身份资格、APEX1／ITSN2反向、QDPR旧源未支持及ARC单比较肽均保留；无统一稳健反向蛋白签名 |
| 空间 | 13标本／5患者；三类配对比较有效患者2／4／3 | 患者级描述、无新P/FDR；RNA派生标签非独立病理真值 |
| GSE122063探索 | 原43在Frontal／Temporal M0有20／17项BH支持；相关供者删除均保持符号 | 持续BH仅13／14；Load Key日期及slide情景为0／0，仍为批次敏感探索；不能按表达修复临床分组 |
| 测量审计 | 固定P1低信号、背景、尺度及供者删除诊断；完整137序列／44蛋白行审计 | 不替换疾病模型；不是新疾病效应、作者实际MSDAP模型复现或功能证明 |

全部阴性、冲突、缺失和不可估计结果按原检验家族保留，不用最小P敏感性、来源换取或测量注释重新筛选81／43。

## 已归档MR及来源限制

正向10工具IVW-MRE为β=−0.06845232870650153，SE=0.017049088100339817，P=5.944224710441849×10⁻⁵；反向57工具β=−0.13533944576480864，SE=0.14133335311078019，P=0.3382694003379726。它们是历史结果的披露，不是本轮重算，也不转为81／43靶点证据。

570个已分析跨方向组合及含全部59个反向暴露工具的590组合均无r²≥0.8，最大仍为0.775069。2项反向结局缺失未补proxy；没有强LD不等于没有局部相关或水平多效性。

Liu补表12／12的R2数值匹配abs-r形式、0／12匹配平方；保留原标签，不宣称作者错误或已获得liability解释度。Liu声明总N24352与Melin分队列24381的版本差异未闭合。0／12实际GBM研究EAF，8／12可追的是1000G EUR参考RAF，不能替代。反向AD为确诊加proxy；**正式liability Steiger仍未运行，无P或PASS**。[逐段来源映射](docs/MR_Paragraph_To_Source_And_Disclosure_Crosswalk.tsv)及`diagnostics/`保留细目。

## 本版实际收录

| 类别 | 实际范围 |
| --- | --- |
| 基因表 | 81／43各1207列；本地1244列原表不变。删除37个可连接身份字段，其余149668个单元逐字相同 |
| 图像 | 43个汇总PNG，含7主图组12组件与31补图；14个供者／Run／个体临床关联图仅本地 |
| 补表 | `supplements/`35个全行原字节汇总文件；211项清单另有2个当前基因表派生映射。其余不是实际payload |
| 复现 | 17条冻结起点范围说明、已归档诊断和文件校验器；没有全部原始输入、17模型实现或环境 |
| 历史实现 | 可读[固定旧提交的实现目录](https://github.com/1209433622cz-maker/MRQI-AD-Glioma/tree/3dd32dbc7d9ce6a6ecbdff8676c5978470556fdc/code/implementation)，仅作历史代码背景；不保证与本版全部17路线兼容 |

旧EDA86／87未发布候选只删14列，漏掉8个最大变化供者及15个测量极值供者字段。本版累计去掉37列，并保存[逐字段去身份日志](docs/IDENTITY_COLUMN_REDACTIONS.tsv)；没有从数值表删基因或不利结果。旧候选不在本版中，Git历史也没有被重写。

17／17不是raw-to-paper自足发布：13条数值／描述、1条原始计数、1条作者汇总提取、2条技术审计。R08历史142／142与当前110／142、32项及385单元差异继续分列；PNG渲染边界不与数值通过混算。复制文件、校验hash和计算复现都不是新增生物学复制。

## 下载与核验

下载本版仓库或固定commit，在任何工作目录运行（Python3，仅标准库）：

```bash
python /path/to/MRQI-AD-Glioma/tools/verify_candidate.py
python /path/to/MRQI-AD-Glioma/tools/verify_submission_alignment.py
```

程序核查精确文件集合、bytes／SHA-256、完整基因身份与1207列、37身份字段排除、43图及35+2补表映射。允许Git本身的`.git`元数据，不将它当发布内容。`PUBLIC_MANIFEST.tsv`不包含自身；`PUBLIC_ASSET_PROVENANCE.tsv`记录逐文件来源。此程序不运行DE、MR或重建原统计模型。

## 版本 权利与引用

版本标识见[PUBLICATION_VERSION.json](PUBLICATION_VERSION.json)，精确commit以所阅Git版本为准；直接父版本为[EDA89固定143文件提交](https://github.com/1209433622cz-maker/MRQI-AD-Glioma/tree/c07775953aa8d395b6b79188d6053427bc0b4a26)；更早66文件版本仍在Git历史。EDA88旧候选与来源审计说明属于其原日期的工作记录，不代表本版当前发布状态。[迁移说明](docs/EDA89_MIGRATION_AND_VERSION_NOTES.md)标出可用性文字修订及保留范围。

项目没有新DOI、正式出版物、MIT／CC0或其他新代码许可证。公开可阅览不等于再许可；应分别引用原研究和固定commit。[权利说明](RIGHTS_AND_SCOPE.md)继续适用。原GEO／PRIDE／PDC数据、第三方全文、原XLSX、参考基因型、个体表、逐人删除记录、受控资料、DUC／签署、账户及凭据不随本版分发。

本版已确定实际55图集合并在本地实装211表附件；其中14图、174表仍仅本地。下一步按目标期刊做连续编号映射和附件要求核对，保持来源和估计对象边界。新的现代GBM对照／分型、GBM精确版本EAF及正式Steiger、因果靶点和真实功能证据仍缺；不能靠公开归档升级结论。本轮不申请新受控资料。

### English summary

This curated EDA91 manuscript-alignment snapshot preserves all81 genes and the original43, negative and conflicting findings, missingness and original multiplicity families. Thirty-seven linkable donor-identity columns are excluded while every retained result cell remains unchanged. It contains43 aggregate images and35 literal whole-table summaries plus2 gene-table derivative mappings;57 figure and211 supplement indexes are not claims of complete payload delivery. Archived MR, source-frequency and sample-version limitations are kept separate from expression evidence. Formal liability-scale Steiger was not run. The snapshot does not establish population protection, same-cell antagonism, a modern GBM-versus-healthy comparison or causal therapeutic targets. The validator checks files rather than replaying all models.
