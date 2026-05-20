# MetaboAnalyst MSEA QEA 输入文件 — 使用指南

## 概述
本目录为 MetaboAnalyst 5.0/6.0 的 **MSEA Quantitative Enrichment Analysis (QEA)**
准备的输入文件. 算法基础是 GlobalTest (Goeman et al. 2004, Bioinformatics),
将每条通路的成员代谢物视为联合预测因子, 在 GLM 框架下检验它们整体是否预测表型.

  Web: https://www.metaboanalyst.ca/
  库: SMPDB (human) + KEGG (human) 双库交叉验证
  模式: QEA (不跑 ORA, 因候选 n=2 功效太弱)
  协变量: age + GA_decimal + year (双哑变量); batch 已在 ComBat 阶段抹掉

## 文件清单

| 文件 | 形状 | 用途 |
|---|---|---|
| **concentration_residualized_80.csv** ★ | 165 × (1+1+67) | **主输入** (80 主轨, 协变量残差化) |
| concentration_residualized_50.csv | 165 × (1+1+73) | 主输入 (50 探索轨, 协变量残差化) |
| concentration_raw_80.csv | 165 × (1+1+67) | 对照 (ComBat 后未残差化) |
| concentration_raw_50.csv | 165 × (1+1+73) | 对照 (ComBat 后未残差化) |
| id_mapping_80.csv | 67 × 6 | 代谢物 → HMDB/KEGG ID 映射 + 缺失诊断 |
| id_mapping_50.csv | 73 × 6 | 同上 |
| samples_metadata.csv | 165 × 7 | 样本元数据 (审计 + MetaboAnalystR 备用) |

浓度表格式: `Sample, Label, Metabolite_1, Metabolite_2, ...`
- Label: `Normal` (n=120) / `Overweight_Obese` (n=45)
- 数值: log2 转换后 (ComBat by injection_batch + 协变量残差化)

## 协变量残差化原理

MetaboAnalyst 5.0 MSEA web 界面不直接支持协变量调整. 我们在外部对每个代谢物
做 OLS 残差化:

```
y_i ~ age + GA_decimal + year_A19 + year_C21    (不含 BMI_group!)
y_i_residualized = residuals + mean(y_i)
```

注意: 协变量回归时 **不放入 BMI_group**, 因此残差矩阵仍保留 BMI 效应,
QEA 跑出的 p 值直接反映 BMI_group 的净效应, 等价于脚本 09 ANCOVA
在通路层面的全样本检验.

## 网页操作流程 (推荐路径)

### Step 1 — 上传数据
1. 访问 https://www.metaboanalyst.ca/, 点击 "Enrichment Analysis"
2. 选 "Concentrations" → 上传 `concentration_residualized_80.csv`
3. Data type: **Continuous** (因为 log2 后是连续值)
4. Sample annotations: 选 "Two-group"

### Step 2 — ID 映射 ("Name Map")
1. MetaboAnalyst 会自动 fuzzy match 代谢物名到 HMDB ID
2. 对未匹配的, 参照 `id_mapping_80.csv` 的 `recommended_input_ID` 列手工补齐
3. 注意: 12-HETrE / 20-HDoHE 等位置异构体 fuzzy match 可能失败, 需手动填 HMDB ID

### Step 3 — 数据预处理 (Data Processing)
**重要**: 因为输入已是 log2 + ComBat + 残差化矩阵, 选:
- Sample normalization: **None**
- Data transformation: **None** (已 log2)
- Data scaling: **None** 或 "Mean centering" (autoscale 会再做一次, 不推荐)

### Step 4 — MSEA QEA
1. 选 "Quantitative Enrichment Analysis"
2. Library 选: **SMPDB (Homo sapiens)** 第一轮, 然后切换 **KEGG (Homo sapiens)** 第二轮
3. Permutations: 2000 (默认) 或 10000 (更稳)
4. 下载 Enrichment Overview Table + bubble plot

### Step 5 — 对照轨 (可选)
重复 Step 1-4, 输入换成 `concentration_raw_80.csv`, 看协变量调整对通路富集的影响

## 论文 Methods 引用模板

> "Pathway enrichment analysis was performed using MetaboAnalyst 5.0
> (Pang et al., Nucleic Acids Res 2021). Concentration matrices (ComBat-corrected
> and covariate-residualized for age, gestational age, and sample year, n=165)
> were submitted to the **Quantitative Enrichment Analysis (QEA)** module, which
> applies the GlobalTest framework (Goeman et al., Bioinformatics 2004) to
> assess whether each pathway's member metabolites jointly predict the BMI
> exposure. Two metabolite set libraries were used in parallel:
> **SMPDB (Homo sapiens)** and **KEGG (Homo sapiens)**.
> Permutation p-values (n=10000) were adjusted by Benjamini-Hochberg FDR.
> Pathway enrichment results were cross-validated against our custom
> enzyme-family GSEA (script 16, see Supplementary)."

## 已知限制 (论文也需明示)

1. **HMDB ID 覆盖率**: 80 主轨 86.6% (58/67), 50 探索轨 86.3% (63/73);
   **关键: 80 主轨两个核心候选 12-HETrE → HMDB0062747, 20-HDoHE → HMDB0060048
   均有 HMDB ID 命中** (虽然 supplier KEGG ID 缺失). KEGG ID 覆盖率仅
   80 主轨 59.7% / 50 探索轨 60.3%, 因此**优先用 HMDB 上传**.
2. **氧化脂质 panel 集中在 hsa00590**: 预期 ~80% 命中 Arachidonic acid metabolism,
   KEGG 区分度天然较差
3. **QEA 不识别异构体方向性**: 位置异构体 (8-HDoHE vs 16-HDoHE) 在 HMDB 可能
   归到同一 ID, 损失部分机制颗粒度 → 这是为什么仍需互校自定义家族 GSEA

## 与已有富集分析的关系

| 层级 | 方法 | 输出 |
|---|---|---|
| L1 自定义 (酶系颗粒度) | 脚本 16 家族 GSEA | results/tables/family_enrichment_ranked.csv |
| L2 KEGG 手工通路 (中等颗粒度) | 脚本 17 | results/tables/kegg_enrichment_{ora,gsea}.csv |
| **L3 MSEA QEA (外部库背书)** | 本目录 + MetaboAnalyst | **待补** |
| L4 三池酶系通路 (机制最细) | 脚本 15g/15h | results/tables/{aa,la}_pathway_ratios_forest.csv |
