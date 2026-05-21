"""19. MetaboAnalyst MSEA QEA 输入文件准备.

为 MetaboAnalyst 5.0/6.0 的 MSEA QEA (Quantitative Enrichment Analysis)
准备三类输入: 浓度矩阵 + ID 映射表 + 样本元数据.

  库选: SMPDB human + KEGG human (双库交叉验证)
  模式: QEA (GlobalTest, 利用全浓度矩阵 + 表型, 不依赖候选筛选)
  协变量: age + GA_decimal + year_A19 + year_C21 (batch 已 ComBat 抹去)

协变量处理 (关键决策, MetaboAnalyst MSEA web 不直接支持协变量):
  对每个代谢物 y_i 做 OLS:
      y_i ~ age + GA_decimal + year_A19 + year_C21   (注: 不含 BMI_group!)
  取 residuals + grand mean, 得到"协变量已洗"的矩阵.
  QEA 跑此矩阵 → p 值直接反映 BMI_group 的净效应, 等价于
  脚本 09 ANCOVA 在通路层面的全样本检验.

  同时输出未残差化的 ComBat 后矩阵作为对照 (raw_*), 方便对比
  "协变量调整对通路富集结果的影响".

数据轨:
  80 主轨 (n=67): 主分析, 与脚本 09/16/17 保持一致
  50 探索轨 (n=73): 探索性, 覆盖度更广 (含 EPA/DHA 几个边缘特征)

输入:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat.xlsx
  data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx
  data/02_preprocessed/sample_alignment_n165.csv

输出 (results/metaboanalyst/):
  concentration_residualized_80.csv   主输入: 165 × 67 协变量残差化矩阵 + Label 列
  concentration_residualized_50.csv   主输入 (探索): 165 × 73
  concentration_raw_80.csv            备用: 未残差化的 ComBat 矩阵 + Label
  concentration_raw_50.csv            备用 (探索)
  id_mapping_80.csv                   67 代谢物 → HMDB/KEGG ID 映射 + 缺失诊断
  id_mapping_50.csv                   73 代谢物 → 映射
  samples_metadata.csv                165 样本元数据 (BMI_group/age/GA/year/batch)
  README.md                           网页操作指南 + 论文 Methods 引用模板
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

ROOT = Path(__file__).resolve().parent.parent
COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
RAW_DIR = ROOT / 'data' / '01_raw'
ALIGN_DIR = ROOT / 'data' / '02_preprocessed'
OUT_DIR = ROOT / 'results' / 'metaboanalyst'
OUT_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

COVARS = ['age', 'GA_decimal', 'year_A19', 'year_C21']

# === MetaboAnalyst 要求 ASCII-only 的代谢物名 ===
# 报错示例: "No special letters (i.e. Latin, Greek) are allowed in feature names!"
# 因此把希腊字母 (α/β/γ/δ/ε/ω, Δ) 替换成 ASCII 拼写; 大小写保持.
GREEK_TO_ASCII = {
    'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta',
    'ε': 'epsilon', 'ω': 'omega', 'Δ': 'Delta', 'Α': 'Alpha',
    'Β': 'Beta', 'Γ': 'Gamma',
}


def asciify_name(name: str) -> str:
    """把代谢物名里的希腊字母替换为 ASCII 拼写, 其它字符不变.

    例: 'α-Linolenic acid' → 'alpha-Linolenic acid'
        'Prostaglandin F2α' → 'Prostaglandin F2alpha'
        '15-Deoxy-δ-12,14-prostaglandin D2' → '15-Deoxy-delta-12,14-prostaglandin D2'
        '11β-13,14-Dihydro-15-keto prostaglandinF2α' →
            '11beta-13,14-Dihydro-15-keto prostaglandinF2alpha'
    """
    if not isinstance(name, str):
        return name
    for greek, ascii_word in GREEK_TO_ASCII.items():
        name = name.replace(greek, ascii_word)
    return name


def load_covariates():
    """加载临床+对齐表, 返回以 omx_id 为索引的协变量表 (与脚本 09 同源).

    含: BMI_group (二分类), age, GA_decimal, year 双哑变量 (B20 ref).
    batch 信息保留但不用于残差化 (已在 ComBat 处理).
    """
    align = pd.read_csv(ALIGN_DIR / 'sample_alignment_n165.csv')
    clin = pd.read_excel(RAW_DIR / '分娩特征分析_2019-2021_主分析队列n165.xlsx')

    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left'))

    # year 派生: A19=2019, B20=2020 (ref), C21=2021
    year_prefix = cov['omx_id'].str[:3]
    cov['year'] = year_prefix
    cov['year_A19'] = (year_prefix == 'A19').astype(float)
    cov['year_C21'] = (year_prefix == 'C21').astype(float)

    return cov.set_index('omx_id')


def load_combat_matrix(track):
    """加载 ComBat 后 log2 矩阵, 返回 (df_meta, df_data: features × samples)."""
    path = COMBAT_DIR / f'ori_n165_filtered{track}_log2_combat.xlsx'
    df = pd.read_excel(path)
    meta = df[[c for c in META_COLS if c in df.columns]].copy()
    # 样本列 = 非 META_COLS & 非 UNIT_COL
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    data = df[['Metabolite Name'] + sample_cols].set_index('Metabolite Name')
    return meta, data


def residualize(data, cov_df):
    """对每个代谢物 (行) 做 OLS y~COVARS 残差化, 加回 grand mean.

    data: features × samples (log2 ComBat 后)
    cov_df: samples × covariates (含 COVARS 列, index = omx_id)
    return: features × samples (残差化后, 量纲与原数据一致)
    """
    # 对齐样本顺序
    common_samples = [s for s in data.columns if s in cov_df.index]
    data = data[common_samples]
    X_raw = cov_df.loc[common_samples, COVARS].values
    X = sm.add_constant(X_raw)

    resid = pd.DataFrame(index=data.index, columns=data.columns, dtype=float)
    for met in data.index:
        y = data.loc[met].values.astype(float)
        # 防御性: 检查是否有 NaN (ComBat 后理论无)
        if np.isnan(y).any():
            raise ValueError(f'NaN found in metabolite "{met}" after ComBat')
        m = sm.OLS(y, X).fit()
        resid.loc[met] = m.resid + y.mean()
    return resid


def build_concentration_table(data_features_x_samples, cov_df, label_col='BMI_group'):
    """构造 MetaboAnalyst 浓度表: rows=samples, col1=Sample, col2=Label, rest=metabolites.

    代谢物名 (列名) 应用 ASCII 化 (asciify_name) 以满足 MetaboAnalyst 要求.
    """
    samples = list(data_features_x_samples.columns)
    # 转置: features × samples → samples × features
    df = data_features_x_samples.T.copy()
    # 代谢物名 ASCII 化 (作为列名)
    df.columns = [asciify_name(c) for c in df.columns]
    table = df.reset_index().rename(columns={'index': 'Sample'})
    labels = cov_df.loc[samples, label_col].values
    # 英文 Label 对 MetaboAnalyst 更稳 (避免中文编码)
    label_map = {'正常': 'Normal', '超重肥胖': 'Overweight_Obese'}
    table.insert(1, 'Label', [label_map.get(l, l) for l in labels])
    return table


def build_id_mapping(meta_df):
    """生成 ID 映射诊断表: 每个代谢物的 HMDB / KEGG ID + ASCII 名 + 缺失标记.

    输出列:
      original_name      — 原始代谢物名 (可能含希腊字母 α/β/δ)
      ascii_name         — ASCII 化代谢物名 (用于浓度表列名 + MetaboAnalyst)
      HMDB ID, KEGG ID
      has_HMDB, has_KEGG
      recommended_input_ID — 推荐 ID 形式: HMDB:xxx / KEGG:xxx / NAME:ascii_name
      asciified          — 是否实际做了 ASCII 替换 (审计用, 便于 reviewer 抽查)
    """
    m = meta_df[['Metabolite Name', 'HMDB ID', 'KEGG ID']].copy()
    m = m.rename(columns={'Metabolite Name': 'original_name'})
    m['ascii_name'] = m['original_name'].apply(asciify_name)
    m['asciified'] = m['original_name'] != m['ascii_name']
    m['has_HMDB'] = m['HMDB ID'].notna() & (m['HMDB ID'].astype(str).str.strip() != '')
    m['has_KEGG'] = m['KEGG ID'].notna() & (m['KEGG ID'].astype(str).str.strip() != '')

    def pick_id(row):
        if row['has_HMDB']:
            return f'HMDB:{row["HMDB ID"]}'
        if row['has_KEGG']:
            return f'KEGG:{row["KEGG ID"]}'
        return f'NAME:{row["ascii_name"]}'
    m['recommended_input_ID'] = m.apply(pick_id, axis=1)

    # 列顺序: 原名 → ASCII 名 → ID → 诊断
    return m[['original_name', 'ascii_name', 'HMDB ID', 'KEGG ID',
              'has_HMDB', 'has_KEGG', 'asciified', 'recommended_input_ID']]


def write_samples_metadata(cov_df):
    """完整 165 样本元数据 (审计 + R API 备用)."""
    # 用 ComBat 后任一矩阵确认样本顺序
    _, data80 = load_combat_matrix(80)
    samples = list(data80.columns)
    md = cov_df.loc[samples, ['BMI_group', 'age', 'GA_decimal', 'year']].copy()
    # batch: 从上机顺序表查 (用于审计, 不入分析)
    try:
        b = pd.read_excel(RAW_DIR / '上机顺序和批次表.xlsx', header=1)
        # 列: injectior, class, sample.name, batch
        b_sub = b[b['class'] == 'Subject'][['sample.name', 'batch']].set_index('sample.name')
        # sample.name 对应 omx_id 还是 标本号? 看脚本 05 ComBat
        # 此处保守: 不强行 merge, 留空让用户确认
        md['injection_batch'] = b_sub.reindex(samples)['batch'].values if all(s in b_sub.index for s in samples) else 'see_combat_audit'
    except Exception as e:
        md['injection_batch'] = f'load_failed: {e}'
    md.index.name = 'Sample'
    label_map = {'正常': 'Normal', '超重肥胖': 'Overweight_Obese'}
    md['Label_EN'] = md['BMI_group'].map(label_map)
    return md.reset_index()


def write_readme():
    text = '''# MetaboAnalyst MSEA QEA 输入文件 — 使用指南

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

## ★ 代谢物名 ASCII 化 (必读)

MetaboAnalyst 拒绝含特殊字符 (拉丁字母、希腊字母 α/β/γ/δ/ω) 的代谢物名,
报错: *"No special letters (i.e. Latin, Greek) are allowed in feature names!"*

因此本目录浓度表的列名 (代谢物名) 全部经过 ASCII 化, 例:
- `α-Linolenic acid` → `alpha-Linolenic acid`
- `Prostaglandin F2α` → `Prostaglandin F2alpha`
- `15-Deoxy-δ-12,14-prostaglandin D2` → `15-Deoxy-delta-12,14-prostaglandin D2`
- `11β-13,14-Dihydro-15-keto prostaglandinF2α` →
  `11beta-13,14-Dihydro-15-keto prostaglandinF2alpha`

`id_mapping_*.csv` 同时保留 **original_name + ascii_name** 两列, 作为本地脚本
(用原名) 和 MetaboAnalyst (用 ASCII 名) 之间的桥梁; `asciified` 列标记是否
实际做了替换 (审计 + reviewer 抽查用).

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
'''
    (OUT_DIR / 'README.md').write_text(text, encoding='utf-8')


def main():
    print('=== 19. MetaboAnalyst MSEA QEA 输入文件准备 ===\n')

    cov = load_covariates()
    print(f'协变量表: {cov.shape[0]} 样本, 含 {list(cov.columns)}')
    print(f'  BMI_group: Normal={int((cov["BMI_group"]=="正常").sum())}, '
          f'Overweight_Obese={int((cov["BMI_group"]=="超重肥胖").sum())}')
    print(f'  year: A19={int((cov["year"]=="A19").sum())}, '
          f'B20={int((cov["year"]=="B20").sum())}, C21={int((cov["year"]=="C21").sum())}')

    for track in [80, 50]:
        print(f'\n--- {track} 轨 ---')
        meta, data = load_combat_matrix(track)
        print(f'  代谢物: {data.shape[0]}, 样本列: {data.shape[1]}')

        # 残差化
        print('  正在做协变量残差化 (OLS y~age+GA+year_A19+year_C21)...')
        data_resid = residualize(data, cov)

        # 浓度表 (主)
        tbl_resid = build_concentration_table(data_resid, cov)
        out_resid = OUT_DIR / f'concentration_residualized_{track}.csv'
        tbl_resid.to_csv(out_resid, index=False, encoding='utf-8-sig')
        print(f'  ✓ {out_resid.relative_to(ROOT)}  ({tbl_resid.shape})')

        # 浓度表 (对照, 未残差化)
        tbl_raw = build_concentration_table(data, cov)
        out_raw = OUT_DIR / f'concentration_raw_{track}.csv'
        tbl_raw.to_csv(out_raw, index=False, encoding='utf-8-sig')
        print(f'  ✓ {out_raw.relative_to(ROOT)}  ({tbl_raw.shape})')

        # ID 映射
        id_map = build_id_mapping(meta)
        out_id = OUT_DIR / f'id_mapping_{track}.csv'
        id_map.to_csv(out_id, index=False, encoding='utf-8-sig')
        n_hmdb = int(id_map['has_HMDB'].sum())
        n_kegg = int(id_map['has_KEGG'].sum())
        n_ascii = int(id_map['asciified'].sum())
        n_total = len(id_map)
        print(f'  ✓ {out_id.relative_to(ROOT)}  HMDB: {n_hmdb}/{n_total} '
              f'({n_hmdb/n_total*100:.1f}%), KEGG: {n_kegg}/{n_total} '
              f'({n_kegg/n_total*100:.1f}%), ASCII-renamed: {n_ascii}')

    # 样本元数据
    md = write_samples_metadata(cov)
    out_md = OUT_DIR / 'samples_metadata.csv'
    md.to_csv(out_md, index=False, encoding='utf-8-sig')
    print(f'\n✓ {out_md.relative_to(ROOT)}  ({md.shape})')

    # README
    write_readme()
    print(f'✓ {(OUT_DIR / "README.md").relative_to(ROOT)}')

    print('\n=== 完成 ===')
    print(f'输出目录: {OUT_DIR.relative_to(ROOT)}')
    print('下一步: 按 README.md 步骤上传 MetaboAnalyst 网页, 优先用 residualized_80.csv')


if __name__ == '__main__':
    main()
