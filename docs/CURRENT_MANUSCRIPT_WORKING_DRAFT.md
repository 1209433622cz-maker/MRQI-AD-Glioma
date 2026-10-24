# 阿尔茨海默病与历史胶质瘤的区室依赖性反向表达证据

RP2.4.1 EDA90 2026-10-24 公开版本与图表引用核对版

Reciprocal expression across tissue compartments in Alzheimer disease and historical glioma

## 中文摘要

评价阿尔茨海默病与历史胶质瘤之间的反向表达证据，区分跨来源方向一致性、区室定位与蛋白测量支持，限定其与疾病因果推断的关系。

整合公开组织、供者级单核、蛋白和肿瘤空间资料，保留完整81／43及原检验家族，分别核对供者、脑区、尺度与测量。将GSE122063既有双通道汇总、供者留一及批次代理结果接入证据链；本轮仅整合，不新增DE或重算FDR。 另整合唯一P1的测量审计及全31供者留一；该诊断不替换原疾病模型。

BA17主分析与单核单基因稳健性未支持统一复制，海马CSRP1提供有限CA1关联，蛋白保留冲突及44来源行测量限制。GSE122063的原43 M0在Frontal／Temporal有20／17项BH支持，均经23次相关供者删除保持符号，持续BH支持仅13／14项；Load Key日期和slide代理下均为0／0。23个多探针原43中8／9个有点方向分歧，20个仅单探针；完整81中AKR1C2不可测。 固定P1低信号压缩经单供者删除仍保留；软件日期与供者／slide／重复混杂。

现有结果支持来源、脑区及历史肿瘤–宿主区室限定的探索性反向表达；GSE122063只能作为批次敏感探索。供者符号稳定不等于批次或测量稳健，未建立统一蛋白签名、人群保护、同细胞拮抗或因果靶点。 高相关和测量留一稳定均不构成尺度等价或新的疾病复制。

## Abstract

To characterize reciprocal expression evidence between Alzheimer disease and historical glioma while separating disease contrasts, tissue compartments and protein measurement support.

Public tissue, donor-level single-nucleus, protein and tumour spatial evidence was integrated while retaining the complete 81-gene and original 43-gene sets and their original testing families. Frozen GSE122063 donor-deletion, probe and batch-proxy results were added without new differential-expression models or FDR calculations in this integration round. A separate fixed-P1 measurement audit and complete donor deletions were integrated without replacing the original disease models.

The BA17 primary family and single-gene single-nucleus sensitivities did not establish uniform replication. Limited CA1 CSRP1 associations and conflicting protein findings remained. In GSE122063, 20/17 of the original 43 genes met BH thresholds in frontal/temporal M0 analyses. All retained their sign in 23 relevant donor deletions, but only 13/14 retained BH support throughout. Load Key date and slide proxies yielded 0/0 BH-supported genes. Of 23 multi-probe genes, 8/9 had a probe with an opposing point estimate; 20 genes had one probe. AKR1C2 remained unmeasured among the 81. Low-signal compression persisted across individual donor deletions; software dates remained confounded with donors, slides and repeats.

Evidence remains bounded by source, brain region and historical tumour-host compartments. GSE122063 provides batch-sensitive exploration, without establishing independent clinical replication, uniform reciprocal protein effects, population protection, same-cell antagonism or causal targets. Sign stability, statistical support and measurement specificity are distinct dimensions. Correlation and measurement stability do not establish scale equivalence or new biological replication.

## Introduction

阿尔茨海默病（AD）与胶质瘤的反向关系提出了一个值得检验的分子问题。既有双向孟德尔随机化研究报告了胶质瘤遗传易感性与较低AD风险的联系[1]，也已有研究结合表达资料及细胞实验提出TNFRSF12A等跨疾病候选[2]。这些结果为研究提供背景，但遗传易感性、人群发病以及患病组织中的表达差异具有不同估计对象，不能相互替代。我们因此将本研究的问题限定为：在明确的来源与比较对象下，是否能识别并解释可追溯的反向表达证据。

组织表达的相反方向未必代表同一种细胞内相反的疾病反应。人类AD的组织及单核研究显示细胞谱系和疾病状态具有不同转录特征[3,4]；GBM单细胞及蛋白组研究则揭示恶性状态、微环境与分子背景的异质性[5,6]。肿瘤核心相对同人宿主脑的差异还包含区室组成与组织背景变化。由此，供者级细胞和空间资料可以帮助定位组织信号，但仅有肿瘤资料时不能估计GBM相对正常组织的细胞疾病效应。

蛋白证据还需区分来源独立性、统计稳健性和测量归属。近期NBB颞叶蛋白研究提供了早发与晚发AD的公开作者汇总和processed资料[7]，使不同脑库家族的方向核对及逐肽审计成为可能。来源全局肽数、比较实际使用肽数和数据库内序列特异性回答不同问题，不应共同换算为一次复制成功。本文在完整保留81／43候选及其阴性、冲突和缺失的基础上，整合区域表达、供者级细胞、蛋白与空间证据，检验反向表达解释能够达到的边界。 海马新增公开来源与逐供者敏感性进一步检验区域表达的可迁移性；平台序列、参考基因与转录本归属则用于区分统计关联和测量特异性，而非替代独立实验。 本次追加GSE122063既有公开探索，重点检验供者方向稳定与批次敏感是否同时存在；它不改变前述验证来源和估计对象。

## Methods

### 研究设计与来源职责

本研究区分AD病例–对照、同人肿瘤–宿主、肿瘤内部状态以及空间区域四类比较。各来源先核对身份、对照、脑区、重复发布及供者层级；同脑库、多平台、同人多区域和同来源敏感性不计为新的独立来源。GSE125583／GSE95587属于同一来源家族；NBB与Emory按已声明脑库家族区分，未作跨脑库基因型去重。所有分析均为对已知公开资料的探索性研究；冻结实施规则不追溯构成预注册。详细分析流程及完整检验家族见补充材料；完整补表索引、实际提供的部分汇总文件、代码和输入校验记录的范围按交付清单区分，不表示公开包包含全部原始数据和模型。EDA17–19属于既有结果之后冻结并执行的公开来源分析；EDA20只作测量与写作整合，不追溯称为盲法预注册。

### 发现检验与顺序候选形成

GSE125583主分析为218AD／70对照，原整数counts采用DESeq2，设计为年龄、性别及诊断；个体RIN、PMI和批次缺失，未声称调整[3,8]。GSE165595经既有全库QC保留30库／15配对患者，含11例历史4级及4例3级；设计为patient加tissue，比较Tumor−同人Normal标签宿主。该宿主标签不等于健康人群对照，亦未统一完成现代IDH-wildtype GBM重分类。两侧经冻结可测过滤与唯一身份桥接形成16,068共同基因。 [9]

R1定义为AD升高／肿瘤相对宿主降低，R2为相反方向。按原双侧P及效应符号取得方向单侧P，以两来源方向P的最大值构成交并检验；在完整32,136方向家族中保留原BH及BY敏感性[10,11]。508候选进入GSE159699的12AD／10Old比较，设计含sex、age、RIN、PMI和diagnosis[12]，原单侧家族仍为508。由此形成81后，GSE13276仅使用5成人历史GBM的肿瘤核心–同人邻近白质，排除全部海绵状血管瘤对照；12数组中的3R／4R在患者内平均。limma配对区室模型的原方向家族为完整81，21项身份缺失保留[13]。81及43均为顺序选择集，不能用其选择来源内的支持率表示独立确认率。 [14]

### 海马公开RNA测序与CA1 CA3芯片验证

GSE278723原counts为19,592行／248列、50供者。50个NX不擅自命名内嗅皮层，一项原始供者编号的位数差异保持原记录，另有一名供者的来源分组冲突不按表达修复；主范围保留191个CA1–CA4样本／48供者。各脑区独立拟合NC、AD、CAA及ADCAA四组，含sex、age、log1p(PMI)，按CPM≥1在≥6样本过滤，TMM及voom后以contrastAsCoef估计四个原比较。主家族43×16=688、81×16=1296及19,592×16=313,472保留；PMI≤48等原敏感性不替代主模型。提供counts的列总量不是FASTQ测序深度，RIN／batch等全部协变量尚不可用。 [15]

GSE29378为63发布数组／33供者，GPL6947已归一化正强度取log2，不重复分位数归一化。原raw检测P及唯一变换指纹闭合数组连接；以冻结HGNC规则和检测P<0.1在≥4/63的探针形成14,507基因。主测量为合格探针log2中位数，最高全样本均值探针仅敏感性，不按疾病P选择。主M0排除一名PMI缺失供者，CA1为31人、CA3为30人，配对CA3−CA1为29人；模型含AD、sex、age、logPMI、bank及chip，配对用chip指标差，并以原平均强度锚定limma方差趋势。完整81有73可测，原43有41可测，三比较家族分别243／129／43,521。 [16]

公开受试者补表在初始芯片效应见到后取得，RIN有61/63值；M4为RIN匹配子集不调整、M5同样本加RIN、M6再限RIN≥5，属于后加敏感性。M5的CA1／CA3／配对单位为31／28／27，配对含RIN均值和差值。不把匹配样本和条件调整混为一项，也不将敏感性q择优替换主q。脑库家族分离不是个体遗传指纹去重。

### 供者 脑库和芯片影响分析

EDA19在该轮效应前固定352个设计，保留完整81／43。GSE29378以M0及后加M5为锚点，在CA1、CA3和同人差分中删除同一供者的全部相关区域；无纳入样本的删除标记no-op，不计独立试验。chip敏感性删除任一区域落在指定chip的整人，簇之间不互斥。分bank仅省略常数bank项，保留其他协变量职责；满秩、残差df≥4且AD／对照各≥3才可估计，不为结果降阶。直接AD×VA交互每锚点另作243／129／43,521校正，未显著不作等效推断。

GSE278723 CA1／CA3按原四组M0删除48名供者，并逐次重做原CPM、TMM及voom。删除系数范围、与基线同号比例及按基线moderated SE标准化的变化仅作描述，不另算最佳删除q；各来源原q原样带入。独立OLS／WLS核验使用R提供的经验贝叶斯方差及自由度，不称第二套eBayes实现。

### 平台序列与公共参考身份审计

EDA20核对完整81相连的全部137条GPL6947平台序列，其中106条为原实际合格探针。2016年GEO注释序列与原平台表逐条比较；针对每条50nt目标及反向互补，按三段seed的pigeonhole保证及全长核对，穷尽检索GENCODE49 ALL转录本和既有GRCh38主组装中0／1／2替换命中。参考分别含533,740条转录本记录、194个主组装序列；官方MD5及本地SHA冻结。GENCODE49 CHR GTF只作当前基因跨度注释，不反推全部ALT或外显子特异性。 [17]

所有命中保留versioned ENST／ENSG和方向。发现不同Ensembl模型同名后，另用GENCODE49官方精确ENST–HGNC表解析，同一HGNC的参考模型不计作新的跨生物基因命中；仅symbol相同不合并。无桥者保留未解析风险。独立BioPython重读全部命中窗口，原81／43及全部44蛋白来源行逐字符串对账。该审计不含indel、个体变异、脑中实际转录本丰度或杂交热力学，不证明全部交叉杂交已排除；仅基因级归属支持不等于全基因转录本覆盖、isoform或蛋白确认。

### BA17与供者级细胞验证

GSE203206使用已发布连续尺度的log2(x+1)，主模型38AD／8对照调整sex、age、RIN，采用HC3及残差df近似t推断；归一化导出对象未完全闭合，不声称重建作者counts流程[18,19]。原43和81家族分别保留，M1–M6及逐供者删除不替代主分析。GSE174367按donor×六谱系聚合整数counts，原DESeq2主设计为Batch、Age、Sex和Diagnosis，RIN及PMI为敏感性[4]。单基因原主家族为43×6=258项、背景486项；R1、R2及AMS43程序另用18项家族。评分为队列／谱系内基因Z值的固定等权组合，不是共同临床量尺，程序支持不下放为组成基因支持。

### 肿瘤状态和空间定位

GSE211376仅使用自身11患者；GSE182109主范围为14名分级一致GBM，2分级冲突及2LGG单列。共同状态在相应患者至少50细胞时比较，先按患者合并区域，再用全源RNA UMI计算log2(CPM+1)；等患者均值峰值、患者一致性及留一分别报告。GSE237183的13标本来自5患者，同患者同区域先聚合；三个固定区域对比合格患者数为2／4／3，报告患者差值及描述范围，不计算新P／FDR。作者RNA派生标签不是独立病理真值，CSRP1／QDPR与标签定义程序的重叠单列；所有tumour-only资料均不承担疾病–正常细胞差异验证。 [20,21,22]

### 蛋白统计和测量核验

既有PXD020296为9AD／10对照，8AsymAD另列，重复Emory供者不双算；CPTAC PDC000204使用11名临床IDH1／2序列检测报告阴性的历史GBM／8同plex GTEx；逐变异及完整IDH2覆盖清单未闭合，不等同于独立WHO2021分型，疾病与采集来源混杂仍在[6]。原完整病例OLS为主、HC3为稳健性，主候选实际可估计BH分母分别24／36；同来源阶段比较的主HC3采用两比较48个可估计项，PMI敏感性及计划86项家族另列，不作为独立复制。 [23]

NBB公开S2含7,388 protein／isoform行，115条metadata中作者排除1条LOAD，保留组为EA11／EOAD25及LA31／LOAD47，不等同每蛋白有效N[7]。主MSqRob、敏感性DEqMS沿用作者原log2FC／P／q和原比较／方法家族，未加入临床协变量，未反推缺失SE或新增候选BH。原单一精确symbol、参考accession base及作者global肽数规则产生的33／22身份资格不改，多行不合并或择优。MS-DAP版本及比较符号依据原冻结源码核对[24]。

对完整DIA-NN processed逐行核验全部44来源entry，以精确Run／文件名连接作者metadata，并将421个未修饰序列在完整作者FASTA中作精确及I/L等价匹配。gene-specific要求全部匹配GN完整且仅属目标gene，entry-specific要求唯一accession；这些注释局限于该FASTA，且与作者专用谱库的推断范围不同。Proteotypic和鉴定q只作技术描述。原作者比较used_peptides另列，不以简单覆盖计数冒充MS-DAP最终筛选；实际模型与fallback路径仍NOT_REPRODUCED。

### 缺失与推断边界

完整81／43的原效应、原P／q、资格、阴性、冲突及缺失均保留；固定家族与原实际可估计家族分别报告，不跨来源、尺度、方法或敏感性合并q。未支持不证明零效应，未匹配／不可估计不当作生物零。不同年龄层一显著一不显著不定义交互，患者描述范围不称置信区间。EDA20–23测量及写作轮不新增DE、P／q、靶点遗传或真实功能结果；纳入的EDA17–19统计为其原冻结输出，也不把旧310遗传阴性逐基因移植为当前81／43的结论。 原43主家族q缺失的背景基因不记作HC3阴性；完整逐基因措辞与原字段对照见补表S24a–h及当前完整81／43表，未将57个跨参考清楚交集转为新检验集。

### 原accession与独立RefSeq版本对照

EDA21以GPL6947原平台表的精确Accession及原GI为历史设计依据，135条有明确版本，2条无版本GenBank记录另按原GI核对，不以2016注释表中的无版本ID代替。当前NCBI人类RNA集15分卷经官方MD5与SHA冻结，含273,576条唯一记录、1,090,237,628nt，NM／NR／XM／XR分别为70,976／24,620／126,866／51,114。完整137条仍按0／1／2替换、双方向、无indel检索；930条命中对应716个精确版本，710条有直接HGNC连接，另6条POC1B-DUSP6 readthrough无冻结HGNC桥，保留独立GeneID及身份不确定。原版本与当前集分开；base查询返回的退役XM不冒称当前记录。全部137另作固定参数、完整50nt的带gap描述性版本比较，只解释差异，不用于准入或修改DE／FDR（S23l–u）。 [25,26]

### 公开双通道探索与多维证据分层

新增GSE122063仅作为公开GEO标签下的批次敏感探索[27]。136条双通道记录精确连接68个物理数组、62个供者–脑区和31名供者（12AD、8VaD、11Control）；同人Frontal／Temporal及Cy3／Cy5不视为独立供者。先在染料内汇总技术重复，再等权平衡染料，原EDA46／47疾病模型沿用作者已归一化log2发布尺度，不二次取log或用raw流程替换。另行固定P1仅用于测量诊断，见下节。M0调整年龄、性别和log1p(PMI)，M1无PMI及最高来源均值探针为既有敏感性。完整81中80可测、原43全可测；AKR1C2无唯一合格注释，保持缺失。 [28]（公开补图S34–S35；本地补图S32–S33未随包；补表S26i–j）

AD两脑区、四项痴呆对比及两项配对脑区交互分别沿用162／86、324／172和162／86的计划81／43检验家族；全来源校正另列。不可测位置用于计划家族校正时P置1，显示仍为NA；BH与BY分开。模型和脑区之间不选最小q，不合并异源q，也不将AD显著而VaD未显著解释为AD特异。EDA48只连接既有输出，未新增DE、P、FDR或来源准入。

供者留一使用31个家族情景，其中每个AD–Control脑区实际有23个相关供者删除；不属于比较的供者情景复用M0。批次敏感性分别使用文件名日期、官方Load Key区块日期与slide的染料平衡贡献代理，并另看去除5729技术重测通道的结果。三种代理不等于真实实验批次，配对设计中无变化的代理不构成独立支持；2项SLIDE VaD–Control残差df不足8仍不可估计。单供者符号稳定、效应衰减与SE膨胀分别解释。（公开补图S36–S37；本地补图S39未随包；补表S26k–m）

158个合格候选探针提供原效应和OLS区间以审计点方向，不计算新的探针P／FDR、不选择有利探针。原43有23个多探针基因与20个单探针基因；相反点方向不直接证明交叉杂交或异构体差异。该GPL16699审计与GSE29378的137条序列／GENCODE／RefSeq审计分别保留，平台内序列不共享不能证明参考基因组或isoform特异性。证据表分开来源职责、原统计、供者、批次、测量及因果边界，不生成综合分数或新的Tier-1命名。（公开补图S38；补表S26k–l）

### 固定P1测量审计及供者技术构成

唯一P1是与原疾病模型分开的测量诊断流程。以68数组／136通道、58,717非控制探针和20,926基因为范围，使用mean前景／median局部背景、normexp saddle与offset50；后续M/A构造不再背景校正或追加offset，联合136通道quantile后重建logG／logR而不再取log。原供者–脑区–染料内平均重复，再等权平衡染料并沿用基因探针median。1255四数组仍隔离；该P1没有代替原EDA46／47作者发布尺度上的疾病模型。（公开补图S41；本地补图S40未随包；补表S28j–m）

EDA52用作者发布矩阵58,717探针跨136通道的精确均值定义10层，分位边界为R type7、相等值归较低层；未按候选、诊断或结果切层。逐探针定义差值、平均差和中心误差，报告Pearson、SD比及RMSE。总MSE分为共同平均差平方、探针间稳定平均差变异、探针内跨通道中心剩余三项；这是代数诊断，不是背景处理的因果贡献或校正矩阵。原厂QC标志只按字面计数，不用于删探针或重加权。（公开补图S42–S43；补表S28n–s）

EDA53对31人精确逐一删除全部FC／TC、染料及技术重复，基因层每次保留60供者–脑区单位，20,926基因均完成31次比较，共648,706行。既有P1值、染料汇总和探针选择不重算。报告所有min–max、最大绝对偏离与并列供者，不设科学稳健准入阈值。完整81／43及158候选探针保留，AKR1C2原缺失不补零。留一范围不是置信区间，也不是疾病效应留一验证。（公开补图S45；本地补图S44未随包；补表S28x–ag）

header明确的软件提取日期、slide和染料用于5／9／2个技术子集。3名供者跨日期，保留全部成员；不将其归为单一湿实验batch。逐探针跨通道与单通道跨探针的统计几何分别报告；所有136通道、62同人脑区、12同染料重复对和68定向数组保留。数值独立核对与原核合成QA只验证实现一致性，不单独增加历史路线覆盖、疾病支持或因果等级。（本地补图S44未随包；补表S28aa–ak）

### 计算复现的起点与覆盖范围

计算复现按原始测序、processed计数或归一化矩阵、已拟合效应表及下游检验分层。R02、R03、R05从冻结整数counts分别重放原4、17、24个DESeq2模型；R06从冻结归一化探针矩阵重放原8个limma模型。每条沿用原模型、样本、过滤、身份、对比、缺失规则和检验家族，先比较全来源再核对完整候选。上述模型重放本身不含原始处理；R01另在EDA79补齐GSE165595原FASTQ至计数，其他来源原始测序及CEL至RMA仍未重建，也不等于在全新安装环境中重建全部软件。

R04从冻结效应表重放16,068基因的32,136项方向IUT及原BH／BY，再核对508→81→43的身份、方向和缺失传递。R02全部4张及R03全部17张重放效应表与R04原输入逐字节一致，构成上游模型与已重放联合检验的精确连接；本轮不重新拟合或重算P／q。原R02的Cook缺失、R03及R05各自收敛规则保持分开；未测、身份缺失、未收敛和方向不支持不合并为零效应。

海马及稳健性R14／R15／R16、R07原7模型与46次留一、R10原8蛋白模型及R12作者汇总提取已完成隔离核对。R12只重放解析和分类，非作者实际MS-DAP／MSqRob模型。R13／R17技术审计、R11空间描述及R08固定条件重放分起点记录，不混计为独立生物学验证（S30c／d）。

R13逐肽／Run／FASTA测量审计与R17序列审计分别按已冻结原核离线隔离重放。R13重新扫描7,154,523行processed及完整FASTA，对24项原输出作字面比较；R17从15份完整RefSeq RNA快照和52个GenBank缓存正文重新执行137查询、原汇总、上下文、描述性拟合及验证，对32项原输出比较。科学函数／类、原参考、样本、准入和判定规则保持不变；合成QA先于真实执行。

R17的1,128条TSV字段规则采用零容差，保留身份、版本、列行顺序、重复记录及NA／空白／0区别。gzip比较完整解压正文与CRC并记录新容器hash；JSON只允许冻结的两个实际UTC及PreSearch代码bytes／SHA绑定变化。R13的24项比较也保留独立字面提取及原测量边界；两项技术重放不代表作者实际MS-DAP／MSqRob模型复现或新增生物学复制。

EDA70在真实运行前绑定输入、代码、运行库、隔离根及两套142合同。从原GSE174367 H5按512列CSC块重聚合，重放原17个DESeq2模型、条件留一及独立R／Python验证。R实测OpenBLAS16、OMP／MKL1，native固定1；原seed、模型、家族和S2_PMI_CC_EX未运行边界不变。

R09固定原13份代码、native单线程、seed20261014及bootstrap2000次，从GSE211376完整counts和GSE182109原矩阵重新聚合全RNA分母。12阶段含独立raw与描述校验；65项保留TSV身份、整数、NA及1e-10缩放容差，JSON键／类型／数组序精确，PNG零像素差，SVG仅归一化日期和一致ID（S31d／e）。

R01沿用原30库、15对、GENCODE27／GRCh38及原STAR与featureCounts配置，在独立输出根从60份冻结FASTQ重新运行比对、计数及技术QC。原183项比较含30库各6项及3项全局输出，保留全基因整数、身份、顺序和缺失规则；STAR仅按原合同排除4个时钟／吞吐字段，QC连续值保留原1e-12缩放容差（S32d／k）。

读者交付采用EDA83只读catalog和17条路线说明，区分kit内部材料与canonical数据根、父包和运行时依赖。文件字节验收不执行模型。R13的24项核心输出比较、21项最终验收及64项原核QA分别记录，不合并分母。EDA84仅核对主文、图例和补表的对应关系，不新增DE或P／FDR。

### 已归档MR的诊断与来源限制

遗传结果作为独立于81／43表达候选的已归档背景证据列示。本轮仅将既有正向10个GBM工具、反向已分析57个AD工具的效应及EDA86工具诊断对齐到全文，不重新拟合MR、重算P值或用MR筛选表达基因。精确rsID交集与1000 Genomes EUR 503名founders参考LD分别审计；强LD排除阈值保持r²≥0.8。已分析集合的570个跨方向组合与含全部59个反向暴露工具的590个组合均保留完整结果，不把后者称为59个已完成结局分析的工具。两项反向结局缺失继续保留，不以proxy补救。（MR与来源附加诊断目录；独立于211项补表索引）

GBM工具的来源审计分别追踪Liu 2024原始补表、其样本声明与Melin 2017原始研究[1,29]。保留补表“R2”原标签及原值，同时以F=(β/SE)²检查其与sqrt(F/(F+N−2))及其平方的数值关系；该检查描述补表数值尺度，不将其自动当作二元疾病liability解释度。样本声明、具体分析版本绑定、参考群体RAF与实际研究EAF分别记账，不互相替代。（MR与来源附加诊断目录）

反向AD暴露GCST90027158包含确诊与proxy病例；其57个精确变异频率及病例／对照分层按原来源审计保留。GBM研究变异EAF和精确版本样本绑定未闭合时，正式二元liability Steiger维持技术未决且未运行，不报告Steiger的P值或PASS，也不以一般人群患病率替换proxy表型来源。该限制不改变既有MR点估计，但约束方向性及人群保护解释。（MR与来源附加诊断目录）

## Results

### 顺序选择形成来源限定候选集

16,068共同基因的32,136方向检验保留508候选，R1为291、R2为217；原BY保留130项。GSE159699在原508家族中支持81项方向，另277项仅方向一致、98项点估计相反、52项不可检验；81中20项保留BY508支持。GSE13276对完整81进一步筛出43，含R1 28项、R2 15项；其余16项仅方向一致、SAP30L相反点估计、21项身份缺失。后续所有来源继续保留完整81／43身份，两个筛选来源内的81/81及43/43不作为无偏复制成功率（主图1、主图2；S02–S05、S14、S24）。

### 海马来源提供有限CA1支持并保留冲突

GSE278723在完整原43中的AD−NC方向一致数为CA1 25、CA2 23、CA3 25、CA4 24，其余分别18／20／18／19项相反。CSRP1 CA1 β=1.4359、95% CI[0.8186,2.0532]，q688=0.01770、q1296=0.03334；CA3同向但q688=0.05299，PMI≤48时CA1同向但q688=0.09581，均不能按阈值择优。ABCA6／ITSN2跨四区反向，ARC若干区上调与旧R2相冲突。CAA与ADCAA仅承担原病理组背景比较，不混称AD复制（主图3B；公开补图S13–S19；本地补图S11–S12未随包；S20）。

GSE29378 M0中完整81在CA1为50同向／23相反／8缺失，CA3为44／29／8；原43为26／15／2及23／18／2。CSRP1 CA1 β=0.4104、CI[0.2251,0.5957]，q129=0.01225、q243=0.02308；CA3及同人区域差未过原候选族，不能称CA1特异。BZW2／PRPF4配对区域差原q129约0.01463／0.01225，但M5中分别约0.301／0.378，且匹配子集已弱化。后加M5 CA3中ACTG1 q243=0.02244、方向反旧R2，M0未支持；ARC与ABCA6在两个海马来源间亦有方向冲突，全部保留（主图3C；公开补图S22–S26；本地补图S20、S21、S27未随包；S21）。

### 供者方向稳定并不证明效应大小或脑库可迁移

352设计中346可估计，6个VA分层因残差自由度不足、其中M5配对还秩亏而不准入；288个计划单供者删除含17个no-op，实际271次改变样本。完整81／43逐次表分别28,512／15,136行。CSRP1 CA1在GSE29378 M0的31个实际删除均正（β范围0.3429–0.4372），M5为0.3234–0.4353，GSE278723 CA1的47次实际删除为1.2057–1.5422，均不翻转。平台系数尺度不同，不合并效应或q。IQGAP3在GSE278723 CA3有6次因原CPM过滤新增的缺测，占位保留，不作生物零（主图4A、本地补图S28、S29未随包；S22）。

GSE29378 CA1 CSRP1整chip供者簇删除仍均正，但M0 β范围0.2452–0.7432，最大变化3.704基线SE，每次删除8–11人。OHSU限定分层仍正，VA不能独立完整估计；两个锚点的直接bank交互均未过预设FDR，不能证明一致或排除混杂。BZW2／PRPF4在单供者及chip簇删除中始终负，故原区域差并非单个极端供者所致，但不能抵消RIN敏感性FDR未过。独立数值复核25,852个可测候选拟合行，明确含no-op（主图4B、补图S30；S22）。

### BA17主分析未提供FDR复制支持

GSE203206主模型的43候选中41项可测，29项与原AD方向一致、12项点估计相反、2项缺失，原BH43均未达阈值；81背景中79项可测，62项方向一致、17项相反，原BH81同样无支持。FOXO4的主q43=0.3425，而相对尺度M6为0.006706，提示尺度依赖，不能替换主分析阴性。41项可测候选中20项在至少一次逐供者删除中改变点估计符号；相反点方向本身也不等于显著反向证据（主图3A；S06）。

### 程序关联与肿瘤定位具有不同稳定性

11AD／7对照的供者级六谱系再分析中，主258项检验有12个基因–谱系BH命中，涉及9基因；BY及加入RIN后均无单基因命中。兴奋性神经元AMS43的主β=1.212，95% CI 0.657–1.768，原BH18 q=0.01183；加入RIN后β=1.047，q=0.02135，程序关联保留。EX的PMI模型仅11AD／4对照，未执行，不记为阴性。这些结果支持特定谱系的程序关联，不能证实同一种细胞在两疾病中的相反反应（主图6A、补图S4；S07）。

两个tumour-only来源的共同三状态比较分别覆盖10和7名患者。43候选中41项可比较，39项具有相同的来源均值峰值，2项不同、2项不可比较；但仅17/43在两来源全部合格患者中一致，32/43在相应全部逐患者删除中保持相同峰值。R1主要定位于少突胶质细胞状态，R2主要定位于恶性状态；区域纳入及删除仍影响个别患者结果。来源均值定位一致不能替代患者一致性或GBM对正常细胞的疾病效应（主图6B；S08）。

### 既有蛋白方向支持未形成双侧稳健确认

原43中，Emory AD的24项可估计蛋白有6项原OLS方向FDR支持，APEX1／ITSN2显著反向；HC3方向支持限于CSRP1、ARC、BZW2。CPTAC的36项可估计蛋白有22项OLS及13项HC3方向支持；ARC／SUN2为OLS显著反向。PSME3、CKAP4、NRIP2在两蛋白来源均有OLS方向支持，却没有候选同时保留两侧HC3支持。相关结果受采集来源、历史诊断及杠杆影响限制（主图5A、补图S2；S09）。

同一Emory来源的AsymAD−Control在base及PMI模型均无原BH支持；AD−AsymAD中CSRP1与ARC的base HC3 q分别为0.005666和0.02510。加入PMI后点方向未翻转，q升至0.1652和0.2980，未维持原BH支持。8名AsymAD中仅2名同时落在三组共同年龄与PMI范围；这一横断面比较不构成纵向进展、PMI因果解释或新独立复制（公开补图S9；本地补图S10未随包；S17）。

### NBB公开汇总提供部分支持和明确冲突

完整81中41基因匹配44蛋白来源行，33项符合原gene身份资格、28项有疾病效应；原43中对应22项及18项，其余单肽、多行和缺匹配状态不删。MSqRob在原43的EOAD比较中有3项期待、3项反向q支持，LOAD有3项期待、1项反向。ARC与QDPR在两年龄层两方法均为期待方向，但ARC每比较只用1条肽。QDPR的NBB效应为+0.7617／+0.3321，旧GBM蛋白方向相符，旧Emory AD却为−0.01870、q=0.90380，未支持。APEX1两年龄层两方法均反向；ITSN2在EOAD两方法及旧Emory反向，PRPF4反向q支持仅MSqRob保留。SUN2／BZW2的年龄层特异阈值结果不构成交互证据（主图5A、主图5B、补图S8；S15–S16）。

### 测量审计与空间描述界定支持范围

完整processed含7,154,523数据行，115/115 Run精确连接作者metadata，无保留样本Run缺失。44来源行中36行具有至少2条I/L gene-specific序列、7行仅1条、1行0条；这是来源行数，不是确认基因数。11行有entry-specific序列，而11条isoform后缀行均无此类观测，不能据此认定isoform缺失。40/44全局数等于作者S2，四行各少1条；176比较槽中140有作者used数，简单双组覆盖数仅21/140相等。ARC两条序列覆盖不同，原四比较used=1保持。完整FASTA特异性、鉴定q、最终比较肽和疾病q仍分开，实际MS-DAP模型未复现（主图5B、公开补图S7；本地补图S6未随包；S18–S19）。 该逐肽／Run／FASTA技术路线已在EDA55按完整原核重放并通过24／24输出比较；技术一致性不补齐作者实际模型日志。

GSE237183空间资料的三个固定区域对比只有2／4／3名合格患者。前两比较均有39/43患者均值方向符合预期，但9/43及25/43存在患者间混合方向；第三比较为26/43均值方向符合预期，37/43混合方向。RNA标签定义与CSRP1／QDPR的重叠限制独立解释；范围是患者点估计范围，未作新P／FDR，不能称为39基因显著复制（补图S3；S10）。

### 全候选序列审计支持归属而不确认异构体

137条50nt序列在GEO原平台表与注释表中完全一致，共得到1,447条转录本及126条主组装命中。官方精确ENST–HGNC连接后，106条原合格探针中94条匹配预期基因且在本检索范围内无跨基因命中，6条仍有不同或未解析基因模型命中，5条无转录本命中，1条仅命中非预期基因。在GENCODE49范围内，完整81按原合格探针组成分为61项参考内归属支持、8项混合／风险、4项缺预期精确命中和8项原不可测；原43相应为33／5／3／2。这是测量注释计数，不是61个疾病复制成功；此处33也不等同原蛋白33／22资格（主图7；S23）。

在GENCODE49范围内，CSRP1唯一合格探针ILMN_1811921匹配同一基因的17／31条参考转录本记录、主组装1个位点，≤2替换范围内未见跨基因转录本命中；它不区分该17条，也不代表31条均在脑中表达。CKAP4／LLGL2／QDPR缺少原合格探针对预期基因的精确转录本匹配，EXOC3L1的原合格探针只命中MATCAP1；这些技术不确定保留，不借参考重映射替换原结果。FLCN、NRIP2、SUN2、SERHL2、NOP58、MYOT及LBH／CLDN9的混合测量风险完整列出。原合格资格、蛋白33／22身份以及44条来源行均未改变（补图S31；S23）。

### 历史设计闭合与当前参考差异并存

137／137序列均在原accession版本或原GI返回记录的声明Probe_Start处精确匹配。原四个未闭合基因的合格探针均位于旧CDS之后：QDPR NM_000320.1、CKAP4 NM_006825.2、LLGL2 XM_946153.1及EXOC3L1 NM_178516.2，原记录均注明后续替换。QDPR与当前NM_000320.3的描述性比对为50个相同碱基及1个query gap，CKAP4与NM_006825.4为50个相同碱基及2个query gap；同分方案分别为4／6，不能唯一定位或称已验证人群indel。LLGL2旧预测记录被NM_004524取代，当前完整集仍无该探针的预期精确匹配。EXOC3L1探针当前仍只精确匹配MATCAP1的NM_001369680.1；不能据此重命名旧候选或重算表达（GENCODE背景见主图7、补图S31；RefSeq细节见补表S23l–u）。

在当前RefSeq范围内，106个原合格探针中96个有预期精确匹配且无其他≤2替换命中，1个同时命中6条未完成HGNC桥接的POC1B-DUSP6 readthrough，9个缺预期精确匹配。完整81对应63项参考内归属清楚、5项混合／风险、5项缺预期精确匹配及8项原不可测；原43为35／3／3／2，与GENCODE49的61／8／4／8及33／5／3／2分别保留。两套参考均清楚的交集为57项，而非将更有利的参考拼成新准入集。RefSeq另提示KIF19、MAGED1、SLC25A48的部分探针不匹配；DUSP6共享readthrough序列的意义未定，不能宣称独立蛋白交叉反应。CSRP1精确匹配4条已知RefSeq记录，无其他≤2替换命中；这与GENCODE49中17／31的范围不同，同样不构成异构体或蛋白确认（GENCODE背景见主图7、补图S31；RefSeq细节见补表S23l–u）。

### 供者方向稳定但批次与测量敏感

GSE122063的M0在原43预定AD家族中，Frontal／Temporal分别20／17项BH支持、5／5项BY支持；同31人两个脑区中的重复支持不算两次独立复制。直接AD–VaD及配对脑区交互未获相应BH支持。所有M0支持项在23次实际相关供者删除中符号稳定（20／20、17／17），但始终保持BH支持的仅13／20、14／17。未删除任何供者以改善主结果。（公开补图S34–S36；补表S26i–l）

原43的BH支持数随批次代理明显变化：M0为20／17，文件名日期为17／18，官方Load Key日期与SLIDE均为0／0，去除5729重测通道为20／17。0／0表示在这些模型及原家族下未达阈值，不表示零效应成立。M0支持子集在Load Key日期下FC／TC中位β比为0.921／0.917、SE比为1.214／1.216；SLIDE下β比为0.484／0.501、SE比为1.379／1.409。该子集仅用于描述衰减与精度，不重定义推断家族或说明批次已被证明为全部原因。（公开补图S37；补表S26k–m）

原43的23个多探针基因中，Frontal／Temporal分别8／9个至少有一个探针点方向与基因M0相反；20个单探针没有探针间复制信息。完整81的AKR1C2保持缺失，不由AKR1C1／3／4替代。所有原探针、阴性、相反和缺失均保留，这些测量限制不因供者留一符号稳定而消失。（公开补图S38；补表S26k–l）

公开Load Key的72数组／144通道／32人中，缺席GEO的一名供者解释了8通道的差额；136／136已发布通道的供者、脑区及染料均精确连接，数组内为同人不同脑区。该桥接未消除论文12AD／9VaD／10Control与GEO12／8／11的临床分组差异。32条通道的文件名日期与Load Key区块日期不一致，后续EDA50已将该区块日期明确为软件FeatureExtractor事件，但真实湿实验事件仍未闭合，因此不将GSE122063升级为正式独立临床验证。（本地补图S33未随包；补表S28a–c）

重校准后的完整81／43按职责保留发现与顺序选择、限定脑区／区室验证、肿瘤内状态、蛋白汇总和GEO标签探索。GSE159699与GSE13276参与81／43形成，BA17主家族未支持，海马有限CSRP1关联及参考依赖不变；蛋白44行、33／22资格、APEX1／ITSN2相反、QDPR旧来源未支持和ARC单比较肽均不变。新增GSE122063只增加批次敏感探索证据，不形成跨尺度一致签名，也未补足GBM健康对照、同细胞疾病效应或因果靶点。

### 低信号压缩与供者技术构成限制

header证据将官方Load Key区块日期对应至FeatureExtractor软件提取事件：68数组均匹配，16数组与文件名日期不同；这没有闭合真实湿实验日期。唯一P1与发布尺度的原43跨62单位Pearson中位数约0.994797，但SD比中位数为0.7281。最低／最高强度层D1／D10的RMSE中位数为3.6775／0.2181 log2，SD比为0.1286／0.9513。高相关因此不能视为幅度或尺度等价。（公开补图S41；本地补图S40未随包；补表S28j–m）

全探针–通道MSE的共同平移平方项、额外探针稳定偏移项和探针内中心剩余项分别占50.2611%、41.5251%和8.2137%。一个常数不能使尺度等价，前两项也不能被称为背景处理的因果解释。31次删除中D1层SD比中位数为0.126934–0.129514，低信号压缩不依赖单个供者；原43逐基因最大绝对变化仍可达Pearson0.041657、SD比0.065222及中心误差RMSE0.081565 log2，故未将所有候选标为稳健。（公开补图S42–S43；补表S28n–s）

5个软件日期子集的D1中心误差中位数为0.8249–1.0073 log2；日期与供者、slide及重复构成混杂，不能认定独立批次贡献。2014-12-03只含三名已有供者的技术重复，具体供者编号仅留本地原表。所有原疾病效应／P／q、44蛋白来源行和33／22资格不变；BA17主分析未支持、蛋白冲突、QDPR旧源未支持及ARC单比较肽限制均未被本轮测量诊断消除。（公开补图S45；本地补图S44未随包；补表S28x–ag）

### 核心计算链条重放保持原结果

发现模型的重放保留R02全部19,108基因×4模型和R03全部21,442基因×17模型；9张与36张科学表分别逐字节一致。R02保留44个来源–模型Cook缺失P，R03有效P无缺失。R04重放保留原32,136方向家族中的12个不可检验位置；21条精确上游连接通过后，原508候选、AD81及限定GBM43的身份与方向均不变。

R05原24模型的17,657来源基因和74张科学表重放一致，12,192候选–模型位置保留10,943可测与1,249缺失；主508仍为81方向支持、277同向未获FDR支持、98冲突和52缺失。R06原8模型的25张科学表一致，完整81的主分析仍为43限定支持、16同向未获FDR支持、1相反及21身份缺失。原43为顺序形成的候选，重放通过不构成新的无偏复制率。

截至EDA79，17条历史路线按限定起点验收17／17：13条数值／描述、1条R01原始计数处理、1项R12字面提取及2项技术审计。R01真实全30库及183项通过；R09原65项通过，AD组成模型未新拟合。R08历史142／142与当前110／142仍分列。全项目raw-to-paper未完成，BA17未支持、蛋白冲突及测量限制不变（S32a–f）。

R13完整原核64项技术检查及24／24输出比较通过；R17原22项验证及32／32严格输出比较通过。R17重新扫描273,576条accession-version记录、1,090,237,628 nt，保留930条命中和716个匹配版本。6个无法精确桥接的版本、9条模糊命中及137／106原查询／准入范围均未改变；这些结果增加技术可追溯性，不改变参考敏感或原测量解释。

EDA60的R11空间描述14／14通过。EDA70的R08从原H5完成17模型、262,092全来源行、完整1458／774及8505条件留一行。历史142／142且文件SHA全部一致；当前控制110／142，15 TSV、17 DDS失败及385个超合同单元格全部保留。缺失、有限数符号与P／FDR 0.05分类无变化；两分支不拼接，原1e-9及PNG零像素容差不变。

R09的45 TSV、12 JSON、5 PNG及3 SVG全部通过；62文件SHA相同，3 SVG仅允许日期／ID差异。完整40,014候选行、两掩码各8,748行及完整81／43保留。缺少原V2冻结清单及空决策目录造成的两次工程停止均留存；补齐后原核、参数及阈值未改，5份部分输出重试前后字节一致。AD组成仍未新拟合（S31d／e）。

R01新30份completion全部通过，未复用历史completion；58,288基因×30库计数、全库QC及过滤结果均通过真实183项比较。终验1821个保护文件无漂移，330条新产物哈希通过；实际历时约43.8小时。新增证据限于原始计数处理的计算可复现性，不新增疾病模型、供者或因果支持（S32d／e／l／n）。

EDA83读者验收覆盖4,872项绑定、3,504个唯一文件，文件差异为0；当前普通WSL与原R03无execute只读guard亦通过。这些结果只确认交付及入口，不重新判定数值结果。R08历史142／142、当前110／142及32项差异继续分列，未据文件验收改判。

### 已归档MR结果及来源审计

已归档的正向10工具IVW-MRE为β=−0.06845232870650153，SE=0.017049088100339817，P=5.944224710441849×10⁻⁵；反向57工具为β=−0.13533944576480864，SE=0.14133335311078019，P=0.3382694003379726。这些是既有分析值的披露，不是本轮新估计，也不构成81／43的靶点证据。已分析集合精确交集为0，570个及扩展590个跨方向组合均无r²≥0.8；最大r²仍为0.775069（rs75061358—rs74504435，距离32976 bp），因此不能称为完全没有局部相关。反向rs12151021、rs12972720缺失结局；精确与强LD排除集合为空，12行既有FE／MRE敏感性结果未改变。（MR与来源附加诊断目录）

12／12个GBM补表“R2”值与sqrt(F/(F+N−2))在预定5×10⁻¹⁰容差内相符，0／12与其平方相符；原XLSX未存公式，不能由此断言作者如何推导或认定作者错误。Liu Table 1及Methods声明6183病例／18169对照（总N=24352），Melin原始分队列表为6191／18190（总N=24381），相差8病例和21对照；具体分析版本连接仍未闭合，未编造QC解释[1,29]。Liu补表未给出12个工具的实际GBM研究EAF；Melin lead表可追到8／12参考RAF，明确来自1000 Genomes EUR，且仅7个属于现用工具，不能替代研究EAF。正式liability Steiger仍未运行。（MR与来源附加诊断目录）

## Discussion

本研究识别了来源、脑区及历史肿瘤–宿主区室限定的反向表达，并在两套公开海马来源中看到CSRP1的有限CA1正向关联。发现层顺序形成81／43并不构成无偏复制；BA17主分析、单核单基因稳健性、蛋白冲突及区域质量依赖仍限制统一外推。单供者方向稳定不抵消批次关联簇下的效应大小敏感性，也不能补足VA脑库精度。完整保留不利结果，是判断当前可支持结论的组成部分。

区室组成是解释组织方向的一种合理工作假说。R1候选在肿瘤少突胶质细胞状态、R2在恶性状态中的定位，与组织对比同时包含宿主和肿瘤成分的情形相容；GBM状态及微环境异质性的既有研究也支持采用这种分层视角[5]。然而，本研究没有在两疾病中估计匹配的同细胞病例–对照效应，也未通过分解或干预证明组成变化导致组织差异。39/41来源均值峰值一致与17/43患者完全一致的区别，尤其说明均值图谱不能替代患者层的复现。

区域与尺度依赖进一步限制了程序的外推。BA17主模型未见原家族FDR支持，FOXO4在相对尺度中的阈值变化提示测量定义会影响结论；它既不能修复主模型阴性，也不能证明脑区之间存在统计交互。供者级程序关联同样不能转化为所有成员基因的疾病支持。上述限制与AD细胞类型和转录状态异质性的原始研究背景相容[3,4]，但本研究尚未将脑区、细胞组成、RIN和其他未测协变量的贡献分离。

海马结果改善了区域表达证据，却未建立脑区特异性或因果保护。CSRP1的CA1效应在单供者删除中稳定，两个来源来自不同脑库家族，但队列年龄、病理、质量信息及测量尺度不同，不能把系数大小平均为共同生物效应。GSE29378的RIN是初始结果之后取得，匹配样本变化、条件调整和精度下降不能混写为单一批次或RIN解释；未显著的bank交互亦不构成等效证明。

蛋白层应理解为逐来源的部分支持，而非统一确认。NBB公开资料扩展了来源家族覆盖[7]，但ARC的比较实际用肽仍为1且旧GBM方向相反；QDPR的NBB与GBM方向组合不消除旧Emory AD未支持；APEX1／ITSN2的反向结果和PRPF4的方法依赖必须同时解释。Emory阶段比较中加入PMI后支持减弱，也不能转写为PMI解释了全部效应。来源间的采集、对照病理背景及估计尺度不同，直接平均效应或累计显著q数会掩盖这些区别。

完整processed及FASTA审计提高了测量透明度，而未新增疾病效应估计。Run连接排除了这一processed表中的保留样本缺失问题；序列归属则显示来源entry、gene和isoform解释并不等同。36条来源行有多条gene-specific序列不是36基因的新复制成功；无entry-specific观测也不是功能isoform不存在。作者专用谱库、全FASTA、鉴定q及比较used肽具有不同范围；计数不一致不能用于调阈值追齐作者结果。缺少实际MS-DAP过滤／fallback日志时，复现边界应继续保留[7,24]。 R13的完整技术重放已闭合本地实现一致性，R17也通过完整参考序列重放；两者均未补充独立供者、功能实验或因果依据。

序列审计使“平台注释单基因”“参考内基因归属”“转录本／蛋白唯一性”有了明确区分。同一HGNC的不同Ensembl模型不应伪增跨基因风险，但其参考位置仍须分别保存；不同或无HGNC连接的gene模型不能仅按相似名称合并。CSRP1的17／31匹配说明它至多提供目标转录本子集的信号，未分离转录本组成变化。61项归属清晰不等于61项疾病关联，也不能把12项混合／不确定或8项缺失悄然删掉后重新校正。

研究的关键局限来自设计和可识别性。发现肿瘤资料包含历史3级与4级，现有资料不能整体重命名为现代IDH-wildtype GBM；同人邻近白质与GTEx也有不同的对照含义。供者较少、协变量重叠不足、RIN／PMI缺失和未闭合的BA17尺度影响精度与外推。空间患者数小且标签由RNA派生，无法独立确认候选功能。来源家族独立性未通过跨脑库基因型验证，已知结果暴露及多来源探索也不具备全新盲法确认的属性。未支持不代表零效应成立，缺失不代表未表达；各原假设家族和不利结果因此完整保留。

历史版本闭合解释了旧探针的设计来源，却不能恢复当代测量专一性。QDPR／CKAP4的短gap差异、LLGL2退役预测转录本与EXOC3L1的当前非预期命中应分别披露，不能统称为错误探针或生物学阴性。RefSeq新增风险而部分GENCODE跨模型风险在RefSeq中未出现，显示参考覆盖和版本的依赖性；两者均非独立实验复制。带gap比对、旧3′UTR归属和HGNC桥只约束解释，不足以确定脑组织中的杂交效率、变异、转录本组成或蛋白水平。完整81／43和全部阴性、冲突、缺失继续保留。 没有主家族q、模型未执行、来源不可测与原FDR未支持应分别表述；对背景基因、蛋白多行、单肽及readthrough的限制不能用归属清楚或方向一致抵消。

这些组织与区室表达发现与既有MR及跨疾病靶点研究处于不同证据层[1,2]。本文现将已归档MR效应、精确及LD诊断和原始来源限制单列披露；本轮没有新拟合MR，也没有完成当前81／43的靶点遗传或真实功能验证。GBM样本版本差异、缺少实际研究EAF、补表R2尺度及反向AD确诊加proxy表型，使正式liability Steiger仍技术未决。跨方向无精确重叠或r²≥0.8不能排除较弱局部LD、水平多效性、选择与生存偏倚，也不能建立现代分型GBM的人群保护效应。旧310基因遗传阴性不能逐项移植至新候选；表达方向不能证明遗传中介。检验人群保护或同细胞机制仍需要相应时间顺序、现代分型、匹配对照与直接功能证据。

计算复现降低实现与传递错误，不弥补来源设计、测量限制或新增供者。R01原始计数处理通过183项，不新增疾病DE。R08固定16线程重现历史142项，当前控制32项差异仍保留。其他原始测序及芯片归一化、B10C神经元组成与NBB实际MS-DAP模型仍有未重放范围。17／17按起点分层，不按阳性强弱排序。

R09已从原counts／RAW重放B10B定位、两CNV掩码及整合。GSE182109主14患者、2分级冲突及2 LGG另列；GSE211376自身11患者，共同三区室7／10。区域和细胞不作患者，CNV不作病理真值；AD组成只读沿用，无新增GBM正常对照或同细胞跨疾病效应（S31d–f）。

GSE122063说明供者影响与技术结构回答不同问题。M0支持项经单人删除仍同向，减少了由某一单人决定符号的疑虑；Load Key日期与slide调整后的原家族0／0则限制了跨批次解释。日期模型中较小β衰减伴随SE增加，slide模型还出现约一半β；小样本、组别与代理结构和额外参数均可能影响估计。现有均值代理不能分离真实批次、病理和组织组成，不能宣布信号全为伪影，也不能依据文件名日期结果更有利而选择它作为结论。

此次证据分层以可回答的问题为单位。原43是顺序选择后的候选，AD组织、历史肿瘤–宿主、CA1／CA3、蛋白及肿瘤内状态的估计对象不同；同向或名义显著不能彼此抵消未支持、相反和缺失。GSE122063保留“批次敏感探索”的限定，不改变既有海马、BA17、蛋白和测量解释，也不以新的公开队列数量替代来源准入与可迁移性检验。

EDA49–53已完成公开临床与日期审计、raw完整性与物理身份核查、唯一P1及供者留一测量诊断；这些任务不再列为未执行。临床分组、真正湿实验事件和作者流程等价仍未闭合，测量支线封版后不继续以变换或删样本恢复疾病阳性。现有独立AD及历史GBM证据需按原职责整合，GBM可解释健康对照和RNA–蛋白异质性仍属科学缺口。

测量审计补充了对尺度依赖的解释，却没有把旧探索提升为独立疾病验证。低信号压缩跨单供者删除保持，可以降低该现象由某一个人决定的疑虑；它不能证明原疾病系数在另一个尺度下成立。单一P1也无法分离前景估计、normexp、offset、quantile及作者未完全复原流程的各自贡献。不得减去估计偏移后制造新矩阵、挑选较有利尺度或据此恢复原家族支持。

本轮将测量技术支线封版：软件日期含义已闭合，真实湿实验事件、Control10／11、AD男性4／3、1255病理与排除身份及作者流程等价仍未闭合。日期子集重叠，不能估计独立日期效应或声称批次已消除。完整证据应回到来源和估计对象：GSE122063保留批次敏感探索，独立脑区、蛋白和肿瘤区室的阴性、冲突与缺失继续限制叙事；没有新增人群保护、同细胞拮抗或因果靶点证据。

## Conclusion

公开资料支持AD与历史胶质瘤之间来源、脑区及肿瘤–宿主区室限定的探索性反向表达和有限CA1关联。GSE122063的原43 M0有20／17项BH支持，Load Key日期和slide代理均为0／0，只能保留批次敏感探索。完整81／43、44蛋白来源行及全部阴性、冲突、缺失和测量限制不变，尚不支持统一蛋白签名、人群保护、同细胞拮抗或因果靶点。 唯一P1测量审计解释了部分尺度限制，未替换疾病模型或改变上述证据等级。

## 数据与代码可追溯性

原始资料来自GEO、PRIDE、PDC与公开参考库，精确编号、冻结文件及来源职责见配套可用性清单。真实公开版本固定为https://github.com/1209433622cz-maker/MRQI-AD-Glioma/tree/c07775953aa8d395b6b79188d6053427bc0b4a26，EDA89已对143个远端文件完成集合、bytes、SHA-256及Git blob逐项验收；此验收不替代模型复现或科学验证。该提交是本文引用的公开证据快照，本地1244列全表与原始输入不等同于该公开快照。

公开快照保留完整81／43各1207列，37个可连接供者身份字段不随包，其余149668个单元与本地1244列原表逐字一致。实际提供43幅原像素汇总图（7组主图12组件及31幅补图）、35项全行原字节汇总表与2项当前基因表派生映射；57图／211补表是完整披露索引，14幅个体／供者／Run关联补图仅本地，不假写为已公开。本文图表引用的公开／本地状态及MR与来源附加诊断14文件目录另列，后者不混入211项补表计数；原始论文、原XLSX及受控资料不随包。

17／17表示各路线冻结起点范围内完成计算验证，公开快照未包含全部raw-to-paper数据、模型和运行时；R08历史与当前环境的32项差异及385个超容差单元仍单列。当前工作稿的EDA90修改仅涉及版本、可用性、图表引用及文字去身份，不新增DE、MR拟合、P／q或因果支持。29条参考文献及各第三方来源条件保留；未新增DOI、代码或稿件开源许可证，不把申请、机构确认或签署材料视为研究数据或作者批准替代。

## References

1. Liu Y, Chen Y, Gao M, et al. Association between glioma and neurodegenerative diseases risk: a two-sample bi-directional Mendelian randomization analysis. Front Neurol. 2024;15:1413015. [原始来源](https://www.frontiersin.org/journals/neurology/articles/10.3389/fneur.2024.1413015/full)

2. Liu T, Pu J, Theil S, et al. Potential role of TNFRSF12A in linking glioblastoma and alzheimer’s disease via shared tumour suppressor pathways. Sci Rep. 2025;15:21535. [原始来源](https://www.nature.com/articles/s41598-025-08000-7)

3. Srinivasan K, Friedman BA, Etxeberria A, et al. Alzheimer’s Patient Microglia Exhibit Enhanced Aging and Unique Transcriptional Activation. Cell Rep. 2020;31:107843. [原始来源](https://doi.org/10.1016/j.celrep.2020.107843)

4. Morabito S, Miyoshi E, Michael N, et al. Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer’s disease. Nat Genet. 2021;53:1143–1155. [原始来源](https://www.nature.com/articles/s41588-021-00894-z)

5. Neftel C, Laffy J, Filbin MG, et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell. 2019;178:835–849.e21. [原始来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/)

6. Wang LB, Karpova A, Gritsenko MA, et al. Proteogenomic and metabolomic characterization of human glioblastoma. Cancer Cell. 2021;39:509–528.e20. [原始来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC8044053/)

7. Thanou E, Ganz A, Pita-Illobre D, et al. Proteomics of post mortem brains in early- and late-onset Alzheimer’s disease: Unraveling differential Aβ effects and potential AD biomarkers. Alzheimers Dement. 2026;22(7):e71662. [原始来源](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.71662)

8. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550. [原始来源](https://doi.org/10.1186/s13059-014-0550-8)

9. Hwang T, Kim S, Chowdhury T, et al. Genome-wide perturbations of Alu expression and Alu-associated post-transcriptional regulations distinguish oligodendroglioma from other gliomas. Communications biology. 2022;5(1):62. [原始来源](https://doi.org/10.1038/s42003-022-03011-w)

10. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B. 1995;57:289–300. [原始来源](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

11. Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency. Ann Stat. 2001;29:1165–1188. [原始来源](https://doi.org/10.1214/aos/1013699998)

12. Nativio R, Lan Y, Donahue G, et al. An integrated multi-omics approach identifies epigenetic alterations associated with Alzheimer’s disease. Nat Genet. 2020;52:1024–1035. [原始来源](https://doi.org/10.1038/s41588-020-0696-0)

13. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47. [原始来源](https://doi.org/10.1093/nar/gkv007)

14. Mangiola A, Saulnier N, De Bonis P, et al. Gene Expression Profile of Glioblastoma Peritumoral Tissue: An Ex Vivo Study. PLoS ONE. 2013;8(3):e57145. [原始来源](https://doi.org/10.1371/journal.pone.0057145)

15. Li C, Chen Y, Luo S, et al. Differential Roles of Astrocytic CSF1 in Alzheimer's Disease and Cerebral Amyloid Angiopathy: Insights from Transcriptomic Analysis. Aging Dis. 2025;16(5):3137–3153. [原始来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC12339160/)

16. Miller JA, Woltjer RL, Goodenbour JM, Horvath S, Geschwind DH. Genes and pathways underlying regional and cell type changes in Alzheimer's disease. Genome Med. 2013;5:48. [原始来源](https://genomemedicine.biomedcentral.com/articles/10.1186/gm452)

17. GENCODE. Human Release 49 GRCh38.p14 and official transcript HGNC metadata. Public reference data. Accessed 18 October 2026. [原始来源](https://www.gencodegenes.org/human/release_49.html)

18. Caldwell AB, Anantharaman BG, Ramachandran S, et al. Transcriptomic profiling of sporadic Alzheimer’s disease patients. Mol Brain. 2022;15:83. [原始来源](https://doi.org/10.1186/s13041-022-00963-2)

19. MacKinnon JG, White H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. J Econom. 1985;29:305–325. [原始来源](https://doi.org/10.1016/0304-4076(85)90158-7)

20. Abdelfattah N, Kumar P, Wang C, et al. Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target. Nature Communications. 2022;13(1):767. [原始来源](https://doi.org/10.1038/s41467-022-28372-y)

21. Ruiz-Moreno C, Salas SM, Samuelsson E, et al. Charting the single-cell and spatial landscape of IDH-wild-type glioblastoma with GBmap. Neuro-Oncology. 2025;27(9):2281-2295. [原始来源](https://doi.org/10.1093/neuonc/noaf113)

22. Greenwald AC, Darnell NG, Hoefflin R, et al. Integrative spatial analysis reveals a multi-layered organization of glioblastoma. Cell. 2024;187(10):2485-2501.e26. [原始来源](https://doi.org/10.1016/j.cell.2024.03.029)

23. Ping L, Kundinger SR, Duong DM, et al. Global quantitative analysis of the human brain proteome and phosphoproteome in Alzheimer’s disease. Scientific Data. 2020;7(1):315. [原始来源](https://doi.org/10.1038/s41597-020-00650-8)

24. Koopmans F. MS-DAP Mass Spectrometry Downstream Analysis Pipeline version 1.0.5. Software. Commit 9d390188c7c91d123975d0fa7befa3bc53a142b2; accessed 24 October 2026. [原始来源](https://github.com/ftwkoopmans/msdap/tree/9d390188c7c91d123975d0fa7befa3bc53a142b2)

25. NCBI RefSeq. Homo sapiens species-specific transcript collection and README. Frozen RNA snapshot acquired 18 October 2026; exact 15-part files and hashes in the source manifest. The species README documents format and is not a versioned RNA release. [原始来源](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/README)

26. NCBI GEO and Nucleotide. GPL6947 original platform accessions and exact-version GenBank records. Public source records accessed 18 October 2026; request response hashes retained. [原始来源](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947)

27. NCBI GEO. GSE122063 public series record and Agilent Microarray Load Key. Frozen public source metadata and normalized series matrix accessed 19 October 2026. [原始来源](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122063)

28. McKay EC, Beck JS, Khoo SK, et al. Peri-Infarct Upregulation of the Oxytocin Receptor in Vascular Dementia. Journal of neuropathology and experimental neurology. 2019;78(5):436-452. [原始来源](https://doi.org/10.1093/jnen/nlz023)

29. Melin BS, Barnholtz-Sloan JS, Wrensch MR, et al. Genome-wide association study of glioma subtypes identifies specific differences in genetic susceptibility to glioblastoma and non-glioblastoma tumors. Nature Genetics. 2017;49(5):789–794. [原始来源](https://doi.org/10.1038/ng.3823)
