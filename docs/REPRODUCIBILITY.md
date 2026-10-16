# 复现层级与依赖

## A. 公开包完整性：可立即运行

Python 3标准库即可：`python tools/verify_public_release.py`。PUBLIC_MANIFEST不列自身，避免自递归；本地报告另外记录外层ZIP和commit。该校验确认资产完整性／固定身份／表结构，不重新拟合DESeq2、limma或蛋白模型。

## B. 实现审阅：已有原代码快照

`code/implementation`收录EDA7实际审计到的统计实现及必要的准备／整合示例。来源路径和哈希见PUBLIC_ASSET_PROVENANCE。原`H:\SCI2\MRQI`路径、历史日期和本地输入依赖保留，以便定位原执行环境；在其他机器重建时须配置root和对应输入。不要仅改变日期后重新执行旧冻结脚本并覆盖旧结果。

| 分析 | 原实现入口 | 依赖／输入层 |
| --- | --- | --- |
| AD counts | rp241_ad_b2e_count_aware.R | 原整数counts、盲QC及年龄/性别设计 |
| 肿瘤counts | rp241_gbm_b1r1_fullgene_discovery.R | 15患者30库featureCounts、固定配对设计 |
| 方向IUT | rp241_ad_gbm_b7_reciprocal_discovery.py | 两发现完整结果与稳定身份；16,068／32,136 |
| AD候选验证 | rp241_ad_gbm_b8_reserved_validation.R | 12AD/10Old、sex/age/RIN/PMI、原508族 |
| 历史GBM限定验证 | rp241_ad_gbm_b9a_gse13276_limited_validation.R | 12数组、患者内平均、原81族 |
| AD谱系程序 | rp241_ad_gbm_b10a_ad_lineage_composition.R | donor×lineage counts、≥50cells、原族与缺失状态 |
| 肿瘤状态 | b10b_count_localization / summarize_state_context | 全源UMI、真实patient/region/author-label关联 |
| 公开蛋白 | b10d_p1_ad_prepare / b10d_p2_protein_direction | 肽身份、GIS/plex与逐基因可估模型 |
| 空间 | b10d_p4_spatial_describe.py | spot作者标签和patient/specimen关联；无新P/FDR |
| BA17 | rp241_eda4_gse203206.py | 发布尺度、七冻结模型、OLS-HC3 |

原计算记录确认R4.1.2、DESeq2 1.34.0、limma3.50.0；量化链STAR2.7.10a、featureCounts2.0.3、GRCh38/GENCODE27、反向链paired fragments。Python部分使用pandas/numpy/scipy/statsmodels等，但部分历史精确版本未完全闭合；不把当前写稿环境伪装为历史分析环境，不提供虚假的完整lockfile。代码片段的依赖脚本并未全部移入公开镜像，完整raw-to-results复跑仍需取得原始数据、设计／身份冻结文件与对应准备代码。

## C. 科学独立确认：尚未完成

包的hash通过、代码可读和原模型可重跑都不是独立来源复制。新GBM对照／独立AD蛋白须另行通过来源、供者家族、模型、权限及真实作者复核。缺失结果保持NA/status；代码中的校正占位1不是实测P。不同源β和Z程序不可合成共同临床量尺。
