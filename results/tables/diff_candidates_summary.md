# 差异代谢物候选筛选结果

**生成时间**: 2026-05-18  
**筛选标准**: `p_limma < 0.05` AND `|log2FC| ≥ log2(1.2) ≈ 0.263`  
**输入**: `results/tables/ancova_main_{80,50}.csv` (ANCOVA v2, 方案 D + 3 协变量)

## 双轨结果对比

| 项 | 80% 主轨 | 50% 探索轨 |
|---|---|---|
| 总特征数 | 67 | 73 |
| 候选数 | **2** | **4** |
| H0 假阳性期望 | 3.35 | 3.65 |
| 实际命中 vs 期望 | 2 vs 3.35 (低于期望 ✓) | 4 vs 3.65 (略高 ⚠) |
| 离群稳健候选 | 2 / 2 | 4 / 4 |
| 方向 (↑ / ↓) | 2 / 0 | 4 / 0 |

## 80% 主轨候选清单 (主分析)

| Metabolite Name | Chinese Name | Class | direction | log2FC | log2FC_CI_lo | log2FC_CI_hi | p_ols_hc3 | p_limma | p_wilcoxon | n_outliers | attenuation_pct | is_robust_to_outlier |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid | 12-羟基-8Z,10E,14Z-二十碳三烯酸 | Fatty Acyls | ↑(超重肥胖>正常) | 0.324 | 0.020 | 0.628 | 0.0369 | 0.0212 | 0.0379 | 1 | -0.082 | True |
| 20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid | 20-羟基-4Z，7Z，10Z，13Z，16Z，18E-二十二碳六烯酸 | Fatty Acyls | ↑(超重肥胖>正常) | 0.319 | 0.039 | 0.600 | 0.0255 | 0.0236 | 0.0313 | 1 | 0.101 | True |

## 50% 探索轨候选清单 (补救 / 敏感性)

| Metabolite Name | Chinese Name | Class | direction | log2FC | log2FC_CI_lo | log2FC_CI_hi | p_ols_hc3 | p_limma | p_wilcoxon | n_outliers | attenuation_pct | is_robust_to_outlier | also_in_80_main |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5,6-DiHydroxy-8Z,11Z,14Z-eicosatrienoic acid | 5，6-双羟基-8Z，11Z，14Z-二十碳三烯酸 | Fatty Acyls | ↑(超重肥胖>正常) | 0.612 | 0.212 | 1.012 | 0.0027 | 0.0190 | 0.0819 | 3 | 0.249 | True | True |
| 12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid | 12-羟基-8Z,10E,14Z-二十碳三烯酸 | Fatty Acyls | ↑(超重肥胖>正常) | 0.324 | 0.020 | 0.628 | 0.0369 | 0.0214 | 0.0379 | 1 | -0.082 | True | True |
| 20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid | 20-羟基-4Z，7Z，10Z，13Z，16Z，18E-二十二碳六烯酸 | Fatty Acyls | ↑(超重肥胖>正常) | 0.319 | 0.039 | 0.600 | 0.0258 | 0.0241 | 0.0313 | 1 | 0.102 | True | True |
| 6,15-Diketo-13,14-dihydro-prostaglandin F1α | 6,15-双酮基-13，14-二氢-前列腺素 F1α | Fatty Acyls | ↑(超重肥胖>正常) | 0.883 | 0.117 | 1.650 | 0.0240 | 0.0474 | 0.0448 | 4 | 0.018 | True | True |

## 双轨重叠分析

- 80 主轨 ∩ 50 探索轨 = **2** 个 (核心一致信号)
- 仅 50 轨独有 = **2** 个 (50% 检出阈值放宽抢救的低丰度信号)

**核心一致候选 (2)**:
  - `12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid` → ↑(超重肥胖>正常) log2FC=0.324, p_limma=0.0212
  - `20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid` → ↑(超重肥胖>正常) log2FC=0.319, p_limma=0.0236

**50 轨独有 (2)**:
  - `5,6-DiHydroxy-8Z,11Z,14Z-eicosatrienoic acid` → ↑(超重肥胖>正常) log2FC=0.612, p_limma=0.0190
  - `6,15-Diketo-13,14-dihydro-prostaglandin F1α` → ↑(超重肥胖>正常) log2FC=0.883, p_limma=0.0474

## Methods 写作要点

- 筛选标准应明示为 **hypothesis-generating tier**, 不作 FDR 控制声明
- 推荐措辞: "Candidate metabolites were defined by `p_limma < 0.05` and `|log2FC| ≥ log2(1.2)`. As no metabolite reached BH-q < 0.05, this threshold serves as a hypothesis-generating tier rather than confirmed differential abundance."
- 每个候选另由 OLS+HC3、Mann-Whitney U、studentized-residual outlier 稳健性互校 (见表)
- 与上游 PERMANOVA (BMI p≈0.40) 与 OPLS-DA (Q² 为负) 的全局阴性一致, 进一步支持 "per-feature evidence is weak but family-level direction is informative"