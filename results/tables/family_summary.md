# 代谢物 → 家族映射汇总

映射表: `data/02_preprocessed/metabolite_family_map.csv`  
生成: `scripts/12_metabolite_family_map.py`

## 7 大家族总览

| 家族 | 总数 | 在 80 主轨 | 80 候选 | 50 候选 |
|---|---|---|---|---|
| **AA-COX** | 23 | 22 | 0 | 1 |
| **AA-CYP/sEH** | 8 | 8 | 0 | 1 |
| **AA-LOX** | 10 | 8 | 1 | 1 |
| **EPA/DHA/DPA-oxylipin** | 15 | 13 | 1 | 1 |
| **Endocannabinoid** | 2 | 2 | 0 | 0 |
| **Free PUFA** | 5 | 5 | 0 | 0 |
| **LA-oxylipin** | 10 | 9 | 0 | 0 |

## 各家族详细成员

### AA-COX  (23 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `12-HHT` | HHT (COX byproduct) | AA | COX/TXA-syn | ✓ | ✓ | - |
| `15d-PGJ2-like` | PG (dehydration) | AA | COX/PGD-syn | ✓ | ✓ | - |
| `δ12-PGJ2` | PG (dehydration) | AA | COX/PGD-syn | ✓ | ✓ | - |
| `Dihomo-PGE2` | PG (homolog) | AdrA | COX | ✓ | ✓ | - |
| `11β-PGE2` | PG (isomer) | AA | COX | ✓ | ✓ | - |
| `11β-dh-keto-PGF2α` | PG (metabolite) | AA | COX/15-PGDH | ✓ | ✓ | - |
| `13,14-dh-15k-PGE2` | PG (metabolite) | AA | COX/15-PGDH | ✓ | ✓ | - |
| `13,14-dh-15k-PGF2α` | PG (metabolite) | AA | COX/15-PGDH | ✓ | ✓ | - |
| `15-keto-PGE2` | PG (metabolite) | AA | COX/15-PGDH | ✓ | ✓ | - |
| `15-keto-PGF2α` | PG (metabolite) | AA | COX/15-PGDH | ✓ | ✓ | - |
| `PGD1` | PG-1 | DGLA | COX/PGD-syn | ✓ | ✓ | - |
| `PGE1` | PG-1 | DGLA | COX/PGE-syn | ✓ | ✓ | - |
| `PGF1α` | PG-1 | DGLA | COX | ✓ | ✓ | - |
| `13,14-dh-15k-PGE1` | PG-1 (metabolite) | DGLA | COX/15-PGDH | - | - | - |
| `13,14-dh-15k-PGF1α` | PG-1 (metabolite) | DGLA | COX/15-PGDH | ✓ | ✓ | - |
| `15-keto-PGE1` | PG-1 (metabolite) | DGLA | COX/15-PGDH | ✓ | ✓ | - |
| `6,15-diketo-dh-PGF1α` | PGI2 metabolite | AA | COX/PGIS/15-PGDH | ✓ | ✓ | ★ |
| `PGA2` | Prostaglandin | AA | COX | ✓ | ✓ | - |
| `PGD2` | Prostaglandin | AA | COX/PGD-syn | ✓ | ✓ | - |
| `PGE2` | Prostaglandin | AA | COX/PGE-syn | ✓ | ✓ | - |
| `PGF2α` | Prostaglandin | AA | COX/PGF-syn | ✓ | ✓ | - |
| `PGJ2` | Prostaglandin | AA | COX/PGD-syn | ✓ | ✓ | - |
| `TXB2` | Thromboxane | AA | COX/TXA-syn | ✓ | ✓ | - |

### AA-CYP/sEH  (8 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `11,12-DiHETrE` | DiHETrE (diol) | AA | sEH | ✓ | ✓ | - |
| `14,15-DiHETrE` | DiHETrE (diol) | AA | sEH | ✓ | ✓ | - |
| `5,6-DiHETrE` | DiHETrE (diol) | AA | sEH | ✓ | ✓ | ★ |
| `8,9-DiHETrE` | DiHETrE (diol) | AA | sEH | ✓ | ✓ | - |
| `11,12-EpETrE` | EpETrE (epoxy) | AA | CYP-Epo | ✓ | ✓ | - |
| `8,9-EpETrE` | EpETrE (epoxy) | AA | CYP-Epo | ✓ | ✓ | - |
| `16-HETE` | HETE-ω | AA | CYP4 | ✓ | ✓ | - |
| `18-HETE` | HETE-ω | AA | CYP4 | ✓ | ✓ | - |

### AA-LOX  (10 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `15-HEDE` | HEDE (20:2 monohydroxy) | 20:2 | 15-LOX | ✓ | ✓ | - |
| `11-HETE` | HETE | AA | 11R-LOX/Auto-ox | ✓ | ✓ | - |
| `12-HETE` | HETE | AA | 12-LOX | ✓ | ✓ | - |
| `15-HETE` | HETE | AA | 15-LOX | ✓ | ✓ | - |
| `8-HETE` | HETE | AA | 8-LOX/Auto-ox | ✓ | ✓ | - |
| `12-HETrE` | HETrE (monohydroxy) | AA→ETrE | 12-LOX via HpETE | ✓ | ★ | ★ |
| `15-HETrE` | HETrE (monohydroxy) | DGLA | 15-LOX | ✓ | ✓ | - |
| `8-HETrE` | HETrE (monohydroxy) | AA→ETrE | 8-LOX/Auto-ox via HpETE | ✓ | ✓ | - |
| `LTB4` | Leukotriene | AA | 5-LOX/LTA4-H | - | - | - |
| `LTE4` | Leukotriene | AA | 5-LOX/LTC4-syn | - | - | - |

### EPA/DHA/DPA-oxylipin  (15 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `19,20-DiHDPA` | DiHDPA (DPA) | DPA | sEH | ✓ | ✓ | - |
| `19,20-EpDPA` | EpDPA (DPA) | DPA | CYP-Epo | ✓ | ✓ | - |
| `10-HDoHE` | HDoHE (DHA) | DHA | LOX | ✓ | ✓ | - |
| `11-HDoHE` | HDoHE (DHA) | DHA | LOX | ✓ | ✓ | - |
| `13-HDoHE` | HDoHE (DHA) | DHA | LOX | ✓ | ✓ | - |
| `14-HDoHE` | HDoHE (DHA) | DHA | LOX | ✓ | ✓ | - |
| `7-HDoHE` | HDoHE (DHA) | DHA | LOX | ✓ | ✓ | - |
| `8-HDoHE` | HDoHE (DHA) | DHA | LOX | - | - | - |
| `16-HDoHE` | HDoHE-ω (DHA) | DHA | CYP-ω | ✓ | ✓ | - |
| `20-HDoHE` | HDoHE-ω (DHA) | DHA | CYP-ω | ✓ | ★ | ★ |
| `11-HEPE` | HEPE (EPA) | EPA | 11R-LOX/Auto-ox | - | - | - |
| `12-HEPE` | HEPE (EPA) | EPA | 12-LOX | ✓ | ✓ | - |
| `15-HpEPE` | HpEPE (EPA) | EPA | 15-LOX | ✓ | ✓ | - |
| `PGF3α` | PG-3 (EPA) | EPA | COX | ✓ | ✓ | - |
| `TXB3` | TX-3 (EPA) | EPA | COX/TXA-syn | ✓ | ✓ | - |

### Endocannabinoid  (2 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `AEA` | N-acylethanolamine | AA | NAPE-PLD | ✓ | ✓ | - |
| `LEA` | N-acylethanolamine | LA | NAPE-PLD | ✓ | ✓ | - |

### Free PUFA  (5 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `AA` | Precursor | AA | - | ✓ | ✓ | - |
| `ALA` | Precursor | ALA | - | ✓ | ✓ | - |
| `CLA` | Precursor | LA | - | ✓ | ✓ | - |
| `LA` | Precursor | LA | - | ✓ | ✓ | - |
| `tLA` | Precursor | tLA | - | ✓ | ✓ | - |

### LA-oxylipin  (10 个)

| short | family_sub | substrate | enzyme | in_80 | cand_80 | cand_50 |
|---|---|---|---|---|---|---|
| `12,13-DiHOME` | DiHOME (diol) | LA | sEH | ✓ | ✓ | - |
| `9,10-DiHOME` | DiHOME (diol) | LA | sEH | ✓ | ✓ | - |
| `12,13-EpOME` | EpOME (epoxy) | LA | CYP-Epo | ✓ | ✓ | - |
| `9,10-EpOME` | EpOME (epoxy) | LA | CYP-Epo | ✓ | ✓ | - |
| `13-HODE` | HODE | LA | 13-LOX/Auto-ox | ✓ | ✓ | - |
| `9-HODE` | HODE | LA | 9-LOX/Auto-ox | ✓ | ✓ | - |
| `13-HpODE` | HpODE | LA | 13-LOX/Auto-ox | ✓ | ✓ | - |
| `9-HpODE` | HpODE | LA | 9-LOX/Auto-ox | - | - | - |
| `13-oxoODE` | oxoODE | LA | 13-PGR/Auto-ox | ✓ | ✓ | - |
| `9-oxoODE` | oxoODE | LA | 13-PGR/Auto-ox | ✓ | ✓ | - |

---
**标记说明**: ✓ = 在 80 主轨内; ★ = 差异代谢物候选 (p_limma<0.05 + |log2FC|≥log2(1.2) + 离群稳健)