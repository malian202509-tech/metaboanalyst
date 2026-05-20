# 代谢物 → 家族映射汇总

映射表: `data/02_preprocessed/metabolite_family_map.csv`  
生成: `scripts/12_metabolite_family_map.py`

## 7 大家族总览

| 家族 | 总数 | 在 80 主轨 | 80 候选 | 50 候选 |
|---|---|---|---|---|
| **AA-COX** | 17 | 17 | 0 | 1 |
| **AA-CYP/sEH** | 8 | 8 | 0 | 1 |
| **AA-LOX** | 7 | 5 | 0 | 0 |
| **DGLA-oxylipin** | 9 | 8 | 1 | 1 |
| **Endocannabinoid** | 2 | 2 | 0 | 0 |
| **Free PUFA** | 5 | 5 | 0 | 0 |
| **LA-oxylipin** | 10 | 9 | 0 | 0 |
| **ω-3 PUFA oxylipins** | 15 | 13 | 1 | 1 |

## 各家族详细成员

### AA-COX  (17 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `12-HHT` | HHT (COX byproduct) | AA | COX/TXA-syn | 5 | none | ✓ | ✓ | - |
| `15d-PGJ2-like` | PG (dehydration) | AA | COX/PGD-syn | 5 | none | ✓ | ✓ | - |
| `δ12-PGJ2` | PG (dehydration) | AA | COX/PGD-syn | 5 | none | ✓ | ✓ | - |
| `Dihomo-PGE2` | PG (homolog) | AdrA | COX | 5 | none | ✓ | ✓ | - |
| `11β-PGE2` | PG (isomer) | AA | COX | 5 | none | ✓ | ✓ | - |
| `11β-dh-keto-PGF2α` | PG (metabolite) | AA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |
| `13,14-dh-15k-PGE2` | PG (metabolite) | AA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |
| `13,14-dh-15k-PGF2α` | PG (metabolite) | AA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |
| `15-keto-PGE2` | PG (metabolite) | AA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |
| `15-keto-PGF2α` | PG (metabolite) | AA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |
| `6,15-diketo-dh-PGF1α` | PGI2 metabolite | AA | COX/PGIS/15-PGDH | 5 | none | ✓ | ✓ | ★ |
| `PGA2` | Prostaglandin | AA | COX | 5 | none | ✓ | ✓ | - |
| `PGD2` | Prostaglandin | AA | COX/PGD-syn | 5 | none | ✓ | ✓ | - |
| `PGE2` | Prostaglandin | AA | COX/PGE-syn | 5 | none | ✓ | ✓ | - |
| `PGF2α` | Prostaglandin | AA | COX/PGF-syn | 5 | none | ✓ | ✓ | - |
| `PGJ2` | Prostaglandin | AA | COX/PGD-syn | 5 | none | ✓ | ✓ | - |
| `TXB2` | Thromboxane | AA | COX/TXA-syn | 5 | none | ✓ | ✓ | - |

### AA-CYP/sEH  (8 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `11,12-DiHETrE` | DiHETrE (diol) | AA | sEH | 5 | none | ✓ | ✓ | - |
| `14,15-DiHETrE` | DiHETrE (diol) | AA | sEH | 5 | none | ✓ | ✓ | - |
| `5,6-DiHETrE` | DiHETrE (diol) | AA | sEH | 5 | none | ✓ | ✓ | ★ |
| `8,9-DiHETrE` | DiHETrE (diol) | AA | sEH | 5 | none | ✓ | ✓ | - |
| `11,12-EpETrE` | EpETrE (epoxy) | AA | CYP-Epo | 5 | none | ✓ | ✓ | - |
| `8,9-EpETrE` | EpETrE (epoxy) | AA | CYP-Epo | 5 | none | ✓ | ✓ | - |
| `16-HETE` | HETE-ω | AA | CYP4 | 4 | minor | ✓ | ✓ | - |
| `18-HETE` | HETE-ω | AA | CYP4 | 4 | minor | ✓ | ✓ | - |

### AA-LOX  (7 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `15-HEDE` | HEDE (20:2 monohydroxy) | 20:2 | 15-LOX | 4 | minor | ✓ | ✓ | - |
| `11-HETE` | HETE | AA | 11R-LOX/Auto-ox | 1 | near-exclusive | ✓ | ✓ | - |
| `12-HETE` | HETE | AA | 12-LOX | 4 | minor | ✓ | ✓ | - |
| `15-HETE` | HETE | AA | 15-LOX | 4 | minor | ✓ | ✓ | - |
| `8-HETE` | HETE | AA | 8-LOX/Auto-ox | 1 | near-exclusive | ✓ | ✓ | - |
| `LTB4` | Leukotriene | AA | 5-LOX/LTA4-H | 5 | none | - | - | - |
| `LTE4` | Leukotriene | AA | 5-LOX/LTC4-syn | 5 | none | - | - | - |

### DGLA-oxylipin  (9 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `12-HETrE` | HETrE (monohydroxy) | DGLA | 12-LOX via HpETE | 5 | none | ✓ | ★ | ★ |
| `15-HETrE` | HETrE (monohydroxy) | DGLA | 15-LOX | 4 | minor | ✓ | ✓ | - |
| `8-HETrE` | HETrE (monohydroxy) | DGLA | 8-LOX/Auto-ox via HpETE | 2 | dominant | ✓ | ✓ | - |
| `PGD1` | PG-1 | DGLA | COX/PGD-syn | 5 | none | ✓ | ✓ | - |
| `PGE1` | PG-1 | DGLA | COX/PGE-syn | 5 | none | ✓ | ✓ | - |
| `PGF1α` | PG-1 | DGLA | COX | 5 | none | ✓ | ✓ | - |
| `13,14-dh-15k-PGE1` | PG-1 (metabolite) | DGLA | COX/15-PGDH | 5 | none | - | - | - |
| `13,14-dh-15k-PGF1α` | PG-1 (metabolite) | DGLA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |
| `15-keto-PGE1` | PG-1 (metabolite) | DGLA | COX/15-PGDH | 5 | none | ✓ | ✓ | - |

### Endocannabinoid  (2 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `AEA` | N-acylethanolamine | AA | NAPE-PLD | 5 | none | ✓ | ✓ | - |
| `LEA` | N-acylethanolamine | LA | NAPE-PLD | 5 | none | ✓ | ✓ | - |

### Free PUFA  (5 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `AA` | Precursor | AA | - | NA | NA | ✓ | ✓ | - |
| `ALA` | Precursor | ALA | - | NA | NA | ✓ | ✓ | - |
| `CLA` | Precursor | LA | - | NA | NA | ✓ | ✓ | - |
| `LA` | Precursor | LA | - | NA | NA | ✓ | ✓ | - |
| `tLA` | Precursor | tLA | - | NA | NA | ✓ | ✓ | - |

### LA-oxylipin  (10 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `12,13-DiHOME` | DiHOME (diol) | LA | sEH | 5 | none | ✓ | ✓ | - |
| `9,10-DiHOME` | DiHOME (diol) | LA | sEH | 5 | none | ✓ | ✓ | - |
| `12,13-EpOME` | EpOME (epoxy) | LA | CYP-Epo | 5 | none | ✓ | ✓ | - |
| `9,10-EpOME` | EpOME (epoxy) | LA | CYP-Epo | 5 | none | ✓ | ✓ | - |
| `13-HODE` | HODE | LA | 13-LOX/Auto-ox | 4 | substantial | ✓ | ✓ | - |
| `9-HODE` | HODE | LA | 9-LOX/Auto-ox | 2 | dominant | ✓ | ✓ | - |
| `13-HpODE` | HpODE | LA | 13-LOX/Auto-ox | 4 | substantial | ✓ | ✓ | - |
| `9-HpODE` | HpODE | LA | 9-LOX/Auto-ox | 2 | dominant | - | - | - |
| `13-oxoODE` | oxoODE | LA | 13-PGR/Auto-ox | 3 | substantial | ✓ | ✓ | - |
| `9-oxoODE` | oxoODE | LA | 13-PGR/Auto-ox | 3 | substantial | ✓ | ✓ | - |

### ω-3 PUFA oxylipins  (15 个)

| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|---|---|
| `19,20-DiHDPA` | DiHDPA (DPA) | DPA | sEH | 5 | none | ✓ | ✓ | - |
| `19,20-EpDPA` | EpDPA (DPA) | DPA | CYP-Epo | 5 | none | ✓ | ✓ | - |
| `14-HDoHE` | HDoHE (DHA) | DHA | 12-LOX (maresin前体) | 4 | minor | ✓ | ✓ | - |
| `10-HDoHE` | HDoHE-non-enz (DHA) | DHA | Non-enzymatic (radical) | 1 | near-exclusive | ✓ | ✓ | - |
| `11-HDoHE` | HDoHE-non-enz (DHA) | DHA | Non-enzymatic (radical) | 1 | near-exclusive | ✓ | ✓ | - |
| `13-HDoHE` | HDoHE-non-enz (DHA) | DHA | Non-enzymatic (radical) | 1 | near-exclusive | ✓ | ✓ | - |
| `16-HDoHE` | HDoHE-non-enz (DHA) | DHA | Non-enzymatic (radical) | 2 | dominant | ✓ | ✓ | - |
| `7-HDoHE` | HDoHE-non-enz (DHA) | DHA | Non-enzymatic (radical) | 1 | near-exclusive | ✓ | ✓ | - |
| `8-HDoHE` | HDoHE-non-enz (DHA) | DHA | Non-enzymatic (radical) | 1 | near-exclusive | - | - | - |
| `20-HDoHE` | HDoHE-ω (DHA) | DHA | CYP4 (ω-3) | 4 | minor | ✓ | ★ | ★ |
| `11-HEPE` | HEPE (EPA) | EPA | 11R-LOX/Auto-ox | 2 | dominant | - | - | - |
| `12-HEPE` | HEPE (EPA) | EPA | 12-LOX | 4 | minor | ✓ | ✓ | - |
| `15-HpEPE` | HpEPE (EPA) | EPA | 15-LOX | 4 | minor | ✓ | ✓ | - |
| `PGF3α` | PG-3 (EPA) | EPA | COX | 5 | none | ✓ | ✓ | - |
| `TXB3` | TX-3 (EPA) | EPA | COX/TXA-syn | 5 | none | ✓ | ✓ | - |

---
**标记说明**:
- ✓ = 在 80 主轨内; ★ = 差异代谢物候选 (p_limma<0.05 + |log2FC|≥log2(1.2) + 离群稳健)
- **ev (enzyme_evidence_level)**:
  - **5** = 结构唯一可达 (PG/TX/epoxide/diol/HHT 等; 自由基不可达终产物)
  - **4** = 酶主导 ≥70%, 少量非酶 (12/15-HETE, 14-/17-/20-HDoHE, 16/18-HETE)
  - **3** = 酶 vs 非酶 ~50:50 (oxoODE, 13-HpODE)
  - **2** = 非酶主导 ~60-70% (**16-HDoHE 2026-05-19 重归此**, 4-HDoHE, 9-HODE, 11-HEPE, 8-HETrE)
  - **1** = 几乎纯非酶 (11/8-HETE, 7/8/10/11/13-HDoHE)
  - **NA** = 前体 PUFA, 不参与氧化分类
- **non_enz**: none / minor / substantial / dominant / near-exclusive / NA

## 16-HDoHE 重归属决策 (2026-05-19)
原 enzyme="CYP-ω", family_sub="HDoHE-ω (DHA)" → 现 enzyme="Non-enzymatic (radical)", family_sub="HDoHE-non-enz (DHA)".
依据:
1. **文献**: VanRollins 2008 J Lipid Res — CYP4-DHA 主产物是 19/20/22-HDHA, **16-HDHA 不在 CYP4 偏好位置 (ω-7)**; Yin & Porter 2005 — DHA 自由基氧化 C15 双烯丙位 H 脱氢 → C16 radical addition → 16-OH + 17E 双键重排 (panel 命名 "16-Hydroxy-...,17E,19Z-" 完全契合该机制).
2. **数据 (实证)**: 单分子 ANCOVA log2FC=+0.173 (p=0.19), 与非酶组 7/10/11/13-HDoHE 中位 (+0.17) 高度同质, 远低于纯 CYP 20-HDoHE (+0.333, p=0.022); 加入"(16+20)/DHA"比值反而稀释纯 CYP 信号 (FC 1.31× vs 单 20 的 1.35×); 加入"(7+10+11+13+16)/DHA"自氧化指数效应几乎不变.