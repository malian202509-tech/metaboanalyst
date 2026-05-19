"""12. 代谢物 → 家族 / 酶系映射 (下游 heatmap/酶比值/GSEA 富集的共同基础).

73 个氧化脂质 (50 探索轨为超集; 80 主轨为 67 子集) 按生物学家族分类:

  7 大家族 (family_main, 论文 heatmap 分块用):
    1. Endocannabinoid    : N-acylethanolamines (LEA, AEA)
    2. Free PUFA          : 前体 (AA, LA, ALA, CLA 等)
    3. AA-COX             : PG / TX / HHT (COX 通路)
    4. AA-LOX             : HETE (5/8/11/12/15-LOX) + LT (5-LOX)
    5. AA-CYP/sEH         : HETrE + EpETrE + DiHETrE + 16/18-HETE (ω羟化)
    6. LA-oxylipin        : HODE/oxoODE/HpODE + EpOME/DiHOME
    7. EPA/DHA/DPA-oxylipin : HEPE/HDoHE/PG3/TX3/EpDPA/DiHDPA (ω3 派生)

子家族 (family_sub) + 上游底物 (substrate) + 主导酶系 (enzyme) 是机制层注释,
供酶活性比值 (§13) 和 GSEA-like 富集 (§14) 直接使用.

分类策略:
  - 手工硬编码映射 (73 个特征, 全部可审计)
  - 优先用 metabolite_name 精确匹配
  - 同时打印未分类警告确保 100% 覆盖

输入:
  data/03_batch_corrected/ori_n165_filtered50_log2_combat.xlsx (73 特征超集)
  data/03_batch_corrected/ori_n165_filtered80_log2_combat.xlsx (67 特征子集, 标 in_main_80)
  results/tables/diff_candidates_50.csv (标 is_candidate_50)
  results/tables/diff_candidates_80.csv (标 is_candidate_80)

输出:
  data/02_preprocessed/metabolite_family_map.csv (73 行 × ~12 列)
  results/tables/family_summary.md (家族成员清单 + 候选标记)
"""
import os
import sys
import tempfile
from pathlib import Path
import pandas as pd

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

ROOT = Path(__file__).resolve().parent.parent
COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
PREP_DIR = ROOT / 'data' / '02_preprocessed'
TABLES = ROOT / 'results' / 'tables'

# (family_main, family_sub, substrate, enzyme, short_label, evidence_level, non_enz)
#
# evidence_level:
#   5 = 结构唯一可达 (PG/TX/epoxide/diol/HHT etc; 自由基机制不可达终产物)
#   4 = 酶主导 ≥70%, 少量非酶 (12/15-HETE, 14-/17-/20-HDoHE, 16-/18-HETE 等)
#   3 = 酶 vs 非酶 ~50:50 (oxoODE, 13-HpODE 等)
#   2 = 非酶主导 ~60-70% (16-HDoHE 重归此, 4-HDoHE, 9-HODE, 11-HEPE 等)
#   1 = 几乎纯非酶 / 自由基 (11-/8-HETE, 7/10/11/13-HDoHE)
#   NA = 前体 PUFA, 不参与
#
# non_enz: "none" / "minor" / "substantial" / "dominant" / "near-exclusive" / "NA"
#
# 重大修订 (2026-05-19, 见 docs §7B.5b):
#   - 16-HDoHE: enzyme "CYP-ω" → "Non-enzymatic (radical)";
#     family_sub "HDoHE-ω (DHA)" → "HDoHE-non-enz (DHA)";
#     evidence_level=2, non_enz="dominant" (实证: log2FC+0.173, 与非酶组同质;
#     VanRollins 2008 J Lipid Res — CYP4-DHA 主产物 19/20/22-HDHA, 16-HDHA 不在 CYP 偏好位置)
FAMILY_MAP = {
    # === 1. Endocannabinoid (NAPE-PLD 酶催化, 结构特异) ===
    'Linoleoyl ethanolamide':                                                              ('Endocannabinoid', 'N-acylethanolamine', 'LA',  'NAPE-PLD',           'LEA',           5, 'none'),
    'Anandamide':                                                                          ('Endocannabinoid', 'N-acylethanolamine', 'AA',  'NAPE-PLD',           'AEA',           5, 'none'),

    # === 2. Free PUFA (precursors, 不参与氧化分类) ===
    'α-Linolenic acid':                                                                    ('Free PUFA', 'Precursor',                'ALA', '-',                   'ALA',          'NA','NA'),
    'Arachidonic acid':                                                                    ('Free PUFA', 'Precursor',                'AA',  '-',                   'AA',           'NA','NA'),
    'Linoelaidic acid':                                                                    ('Free PUFA', 'Precursor',                'tLA', '-',                   'tLA',          'NA','NA'),
    'Linoleic acid':                                                                       ('Free PUFA', 'Precursor',                'LA',  '-',                   'LA',           'NA','NA'),
    'Conjugated linoleic acids':                                                           ('Free PUFA', 'Precursor',                'LA',  '-',                   'CLA',          'NA','NA'),

    # === 3. AA-COX (PG/TX/HHT 均结构特异, 自由基不可达) ===
    'Thromboxane B2':                                                                      ('AA-COX', 'Thromboxane',                 'AA',  'COX/TXA-syn',        'TXB2',          5, 'none'),
    '11β-Prostaglandin E2':                                                                ('AA-COX', 'PG (isomer)',                 'AA',  'COX',                '11β-PGE2',      5, 'none'),
    '15-Deoxy-δ-12,14-prostaglandin D2':                                                   ('AA-COX', 'PG (dehydration)',            'AA',  'COX/PGD-syn',        '15d-PGJ2-like', 5, 'none'),
    '11β-13,14-Dihydro-15-keto prostaglandinF2α':                                          ('AA-COX', 'PG (metabolite)',             'AA',  'COX/15-PGDH',        '11β-dh-keto-PGF2α', 5, 'none'),
    '1α,1β-Dihomo prostaglandin E2':                                                       ('AA-COX', 'PG (homolog)',                'AdrA','COX',                 'Dihomo-PGE2',   5, 'none'),
    'Prostaglandin A2':                                                                    ('AA-COX', 'Prostaglandin',               'AA',  'COX',                'PGA2',          5, 'none'),
    'Prostaglandin D2':                                                                    ('AA-COX', 'Prostaglandin',               'AA',  'COX/PGD-syn',        'PGD2',          5, 'none'),
    'Prostaglandin E2':                                                                    ('AA-COX', 'Prostaglandin',               'AA',  'COX/PGE-syn',        'PGE2',          5, 'none'),
    'Prostaglandin F2α':                                                                   ('AA-COX', 'Prostaglandin',               'AA',  'COX/PGF-syn',        'PGF2α',         5, 'none'),
    'Prostaglandin J2':                                                                    ('AA-COX', 'Prostaglandin',               'AA',  'COX/PGD-syn',        'PGJ2',          5, 'none'),
    'δ-12-Prostaglandin J2':                                                               ('AA-COX', 'PG (dehydration)',            'AA',  'COX/PGD-syn',        'δ12-PGJ2',      5, 'none'),
    '15-Keto prostaglandin F2α':                                                           ('AA-COX', 'PG (metabolite)',             'AA',  'COX/15-PGDH',        '15-keto-PGF2α', 5, 'none'),
    '6,15-Diketo-13,14-dihydro-prostaglandin F1α':                                         ('AA-COX', 'PGI2 metabolite',             'AA',  'COX/PGIS/15-PGDH',   '6,15-diketo-dh-PGF1α', 5, 'none'),
    '13,14-Dihydro-15-keto prostaglandin E2':                                              ('AA-COX', 'PG (metabolite)',             'AA',  'COX/15-PGDH',        '13,14-dh-15k-PGE2', 5, 'none'),
    '15-Keto prostaglandin E1':                                                            ('AA-COX', 'PG-1 (metabolite)',           'DGLA','COX/15-PGDH',        '15-keto-PGE1',  5, 'none'),
    '13,14-Dihydro-15-keto Prostaglandin E1':                                              ('AA-COX', 'PG-1 (metabolite)',           'DGLA','COX/15-PGDH',        '13,14-dh-15k-PGE1', 5, 'none'),
    'Prostaglandin F1α':                                                                   ('AA-COX', 'PG-1',                        'DGLA','COX',                'PGF1α',         5, 'none'),
    'Prostaglandin D1':                                                                    ('AA-COX', 'PG-1',                        'DGLA','COX/PGD-syn',        'PGD1',          5, 'none'),
    '13,14-Dihydro-15-keto-pgf2α':                                                         ('AA-COX', 'PG (metabolite)',             'AA',  'COX/15-PGDH',        '13,14-dh-15k-PGF2α', 5, 'none'),
    '15-Keto prostaglandin E2':                                                            ('AA-COX', 'PG (metabolite)',             'AA',  'COX/15-PGDH',        '15-keto-PGE2',  5, 'none'),
    '13,14-Dihydro-15-keto Prostaglandin F1α':                                             ('AA-COX', 'PG-1 (metabolite)',           'DGLA','COX/15-PGDH',        '13,14-dh-15k-PGF1α', 5, 'none'),
    'Prostaglandin E1':                                                                    ('AA-COX', 'PG-1',                        'DGLA','COX/PGE-syn',        'PGE1',          5, 'none'),
    '12S-Hydroxy-5Z,8E,10E-heptadecatrienoic acid':                                        ('AA-COX', 'HHT (COX byproduct)',         'AA',  'COX/TXA-syn',        '12-HHT',        5, 'none'),

    # === 4. AA-LOX (HETE + HETrE + LT) ===
    '11-Hydroxy-5Z,8Z,11E,14Z-eicosatetraenoic acid':                                      ('AA-LOX', 'HETE',                        'AA',  '11R-LOX/Auto-ox',    '11-HETE',       1, 'near-exclusive'),
    '15-Hydroxy-5Z,8Z,11Z,13E-eicosatetraenoic acid':                                      ('AA-LOX', 'HETE',                        'AA',  '15-LOX',             '15-HETE',       4, 'minor'),
    '8-Hydroxy-5Z,9E,11Z,14Z-eicosatetraenoic acid':                                       ('AA-LOX', 'HETE',                        'AA',  '8-LOX/Auto-ox',      '8-HETE',        1, 'near-exclusive'),
    '12-Hydroxy-5Z,8Z,10E,14Z-eicosatetraenoic acid':                                      ('AA-LOX', 'HETE',                        'AA',  '12-LOX',             '12-HETE',       4, 'minor'),
    'Leukotriene E4':                                                                      ('AA-LOX', 'Leukotriene',                 'AA',  '5-LOX/LTC4-syn',     'LTE4',          5, 'none'),
    'Leukotriene B4':                                                                      ('AA-LOX', 'Leukotriene',                 'AA',  '5-LOX/LTA4-H',       'LTB4',          5, 'none'),
    '8-Hydroxy-9E,11Z,14Z-eicosatrienoic acid':                                            ('AA-LOX', 'HETrE (monohydroxy)',     'AA→ETrE','8-LOX/Auto-ox via HpETE', '8-HETrE',  2, 'dominant'),
    '12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid':                                           ('AA-LOX', 'HETrE (monohydroxy)',     'AA→ETrE','12-LOX via HpETE',        '12-HETrE', 5, 'none'),
    '15-Hydroxy-8Z,11Z,13E-eicosatrienoic acid':                                           ('AA-LOX', 'HETrE (monohydroxy)',     'DGLA',   '15-LOX',                  '15-HETrE', 4, 'minor'),
    '15-Hydroxy-11Z,13E-eicosadienoic acid':                                               ('AA-LOX', 'HEDE (20:2 monohydroxy)', '20:2',   '15-LOX',                  '15-HEDE',  4, 'minor'),

    # === 5. AA-CYP/sEH (HETE-ω / EpETrE / DiHETrE) ===
    '16-Hydroxy-5Z,8Z,11Z,14Z-eicosatetraenoic acid':                                      ('AA-CYP/sEH', 'HETE-ω',                  'AA',  'CYP4',                '16-HETE',       4, 'minor'),
    '18-Hydroxy-5Z,8Z,11Z,14Z-eicosatetraenoic acid':                                      ('AA-CYP/sEH', 'HETE-ω',                  'AA',  'CYP4',                '18-HETE',       4, 'minor'),
    '8,9-Dpoxy-5Z,11Z,14Z-eicosatrienoic acid':                                            ('AA-CYP/sEH', 'EpETrE (epoxy)',          'AA',  'CYP-Epo',             '8,9-EpETrE',    5, 'none'),
    '11,12-Epoxy-5Z,8Z,14Z-eicosatrienoic acid':                                           ('AA-CYP/sEH', 'EpETrE (epoxy)',          'AA',  'CYP-Epo',             '11,12-EpETrE',  5, 'none'),
    '5,6-DiHydroxy-8Z,11Z,14Z-eicosatrienoic acid':                                        ('AA-CYP/sEH', 'DiHETrE (diol)',          'AA',  'sEH',                 '5,6-DiHETrE',   5, 'none'),
    '8,9-DiHydroxy-5Z,11Z,14Z-eicosatrienoic acid':                                        ('AA-CYP/sEH', 'DiHETrE (diol)',          'AA',  'sEH',                 '8,9-DiHETrE',   5, 'none'),
    '11,12-DiHydroxy-5Z,8Z,14Z-eicosatrienoic acid':                                       ('AA-CYP/sEH', 'DiHETrE (diol)',          'AA',  'sEH',                 '11,12-DiHETrE', 5, 'none'),
    '14,15-DiHydroxy-5Z,8Z,11Z-eicosatrienoic acid':                                       ('AA-CYP/sEH', 'DiHETrE (diol)',          'AA',  'sEH',                 '14,15-DiHETrE', 5, 'none'),

    # === 6. LA-oxylipin (HODE/HpODE/oxoODE + EpOME/DiHOME) ===
    '9-Hydroxy-10E,12Z-octadecadienoic acid':                                              ('LA-oxylipin', 'HODE',                   'LA',  '9-LOX/Auto-ox',       '9-HODE',        2, 'dominant'),
    '13-Hydroxy-9Z,11E-octadecadienoic acid':                                              ('LA-oxylipin', 'HODE',                   'LA',  '13-LOX/Auto-ox',      '13-HODE',       4, 'substantial'),
    '9-Oxo-10E,12Z-octadecadienoic acid':                                                  ('LA-oxylipin', 'oxoODE',                 'LA',  '13-PGR/Auto-ox',      '9-oxoODE',      3, 'substantial'),
    '13-Oxo-9Z,11E-octadecadienoicacid':                                                   ('LA-oxylipin', 'oxoODE',                 'LA',  '13-PGR/Auto-ox',      '13-oxoODE',     3, 'substantial'),
    '9-Hydroperoxy-10E,12E-octadecadienoic acid':                                          ('LA-oxylipin', 'HpODE',                  'LA',  '9-LOX/Auto-ox',       '9-HpODE',       2, 'dominant'),
    '13-Hydroperoxy-9Z,11E-octadecadienoic acid':                                          ('LA-oxylipin', 'HpODE',                  'LA',  '13-LOX/Auto-ox',      '13-HpODE',      4, 'substantial'),
    '9,10-Epoxy-12Z-octadecenoic acid':                                                    ('LA-oxylipin', 'EpOME (epoxy)',          'LA',  'CYP-Epo',             '9,10-EpOME',    5, 'none'),
    '12,13-Epoxy-9Z-octadecenoic acid':                                                    ('LA-oxylipin', 'EpOME (epoxy)',          'LA',  'CYP-Epo',             '12,13-EpOME',   5, 'none'),
    '9,10-DiHydroxy-12Z-octadecenoic acid':                                                ('LA-oxylipin', 'DiHOME (diol)',          'LA',  'sEH',                 '9,10-DiHOME',   5, 'none'),
    '12,13 -DiHydroxy-9Z-octadecenoic acid':                                               ('LA-oxylipin', 'DiHOME (diol)',          'LA',  'sEH',                 '12,13-DiHOME',  5, 'none'),

    # === 7. EPA/DHA/DPA-oxylipin (ω3 系列) ===
    # EPA
    '12-Hydroxy-5,8,10,14,17-eicosapentaenoic acid':                                       ('EPA/DHA/DPA-oxylipin', 'HEPE (EPA)',     'EPA', '12-LOX',              '12-HEPE',       4, 'minor'),
    '15-Hydroperoxy-5,8,11,14,17-eicosapentaenoic acid':                                   ('EPA/DHA/DPA-oxylipin', 'HpEPE (EPA)',    'EPA', '15-LOX',              '15-HpEPE',      4, 'minor'),
    '11-Hydroxy- 5Z,8Z,12E,14Z,17Z-eicosapentaenoic acid':                                 ('EPA/DHA/DPA-oxylipin', 'HEPE (EPA)',     'EPA', '11R-LOX/Auto-ox',     '11-HEPE',       2, 'dominant'),
    'Prostaglandin F3α':                                                                   ('EPA/DHA/DPA-oxylipin', 'PG-3 (EPA)',     'EPA', 'COX',                 'PGF3α',         5, 'none'),
    'Thromboxane B3':                                                                      ('EPA/DHA/DPA-oxylipin', 'TX-3 (EPA)',     'EPA', 'COX/TXA-syn',         'TXB3',          5, 'none'),
    # DHA — 注意 16-HDoHE 2026-05-19 重归为非酶 (见顶部注释)
    '14-Hydroxy-4Z,7Z,10Z,12E,16Z,19Z-docosahexaenoic acid':                               ('EPA/DHA/DPA-oxylipin', 'HDoHE (DHA)',    'DHA', '12-LOX (maresin前体)','14-HDoHE',     4, 'minor'),
    '8-Hydroxy-4Z,6E,10Z,13Z,16Z,19Z-docosahexaenoic acid':                                ('EPA/DHA/DPA-oxylipin', 'HDoHE-non-enz (DHA)', 'DHA', 'Non-enzymatic (radical)', '8-HDoHE', 1, 'near-exclusive'),
    '10-Hydroxy-4Z,7Z,11E,13Z,16Z,19Z-docosahexaenoic acid':                               ('EPA/DHA/DPA-oxylipin', 'HDoHE-non-enz (DHA)', 'DHA', 'Non-enzymatic (radical)', '10-HDoHE', 1, 'near-exclusive'),
    '11-Hydroxy-4Z,7Z,9E,13Z,16Z,19Z-docosahexaenoic acid':                                ('EPA/DHA/DPA-oxylipin', 'HDoHE-non-enz (DHA)', 'DHA', 'Non-enzymatic (radical)', '11-HDoHE', 1, 'near-exclusive'),
    '13-Hydroxy-4Z,7Z,10Z,14E,16Z,19Z-docosahexaenoic acid':                               ('EPA/DHA/DPA-oxylipin', 'HDoHE-non-enz (DHA)', 'DHA', 'Non-enzymatic (radical)', '13-HDoHE', 1, 'near-exclusive'),
    '7-Hydroxy-4Z,8E,10Z,13Z,16Z,19Z-docosahexaenoic acid':                                ('EPA/DHA/DPA-oxylipin', 'HDoHE-non-enz (DHA)', 'DHA', 'Non-enzymatic (radical)', '7-HDoHE',  1, 'near-exclusive'),
    '20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid':                               ('EPA/DHA/DPA-oxylipin', 'HDoHE-ω (DHA)',  'DHA', 'CYP4 (ω-3)',          '20-HDoHE',      4, 'minor'),
    '16-Hydroxy-4Z,7Z,10Z,13Z,17E,19Z-docosahexaenoic acid':                               ('EPA/DHA/DPA-oxylipin', 'HDoHE-non-enz (DHA)', 'DHA', 'Non-enzymatic (radical)', '16-HDoHE', 2, 'dominant'),
    # DPA
    '19(20)-Epoxy-4Z,7Z,10Z,13Z,16Z-docosapentaenoic acid':                                ('EPA/DHA/DPA-oxylipin', 'EpDPA (DPA)',    'DPA', 'CYP-Epo',             '19,20-EpDPA',   5, 'none'),
    '19,20-DiHydroxy-4Z,7Z,10Z,13Z,16Z-docosapentaenoic acid':                             ('EPA/DHA/DPA-oxylipin', 'DiHDPA (DPA)',   'DPA', 'sEH',                 '19,20-DiHDPA',  5, 'none'),
}


def atomic_write_csv(df, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=target.stem + '_', suffix='.csv.tmp', dir=str(target.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False, encoding='utf-8-sig')
        os.replace(tmp, target)
    except PermissionError as e:
        if os.path.exists(tmp): os.unlink(tmp)
        raise PermissionError(f'无法写入 {target.name}: 文件可能正被 Excel/WPS 打开.') from e
    except Exception:
        if os.path.exists(tmp):
            try: os.unlink(tmp)
            except OSError: pass
        raise


def main():
    print('=== 12. 代谢物→家族 映射 ===\n')

    df50 = pd.read_excel(COMBAT_DIR / 'ori_n165_filtered50_log2_combat.xlsx')
    df80 = pd.read_excel(COMBAT_DIR / 'ori_n165_filtered80_log2_combat.xlsx')

    feat50 = set(df50['Metabolite Name'])
    feat80 = set(df80['Metabolite Name'])
    print(f'  50 轨特征: {len(feat50)}')
    print(f'  80 轨特征: {len(feat80)}  (子集)')

    # 候选标记
    cand80 = pd.read_csv(TABLES / 'diff_candidates_80.csv', encoding='utf-8-sig')
    cand50 = pd.read_csv(TABLES / 'diff_candidates_50.csv', encoding='utf-8-sig')
    cand80_set = set(cand80['Metabolite Name'])
    cand50_set = set(cand50['Metabolite Name'])

    rows = []
    unmapped = []
    for _, r in df50.iterrows():
        name = r['Metabolite Name']
        if name not in FAMILY_MAP:
            unmapped.append(name)
            continue
        fmain, fsub, sub, enz, short, ev, ne = FAMILY_MAP[name]
        rows.append({
            'Metabolite Name':            name,
            'Chinese Name':               r['Chinese Name'] if pd.notna(r['Chinese Name']) else '',
            'short_label':                short,
            'family_main':                fmain,
            'family_sub':                 fsub,
            'substrate':                  sub,
            'enzyme':                     enz,
            'enzyme_evidence_level':      ev,
            'non_enzymatic_contribution': ne,
            'KEGG ID':                    r['KEGG ID'] if pd.notna(r['KEGG ID']) else '',
            'HMDB ID':                    r['HMDB ID'] if pd.notna(r['HMDB ID']) else '',
            'in_main_80':                 name in feat80,
            'in_explor_50':               True,
            'is_candidate_80':            name in cand80_set,
            'is_candidate_50':            name in cand50_set,
        })

    if unmapped:
        print(f'\n⚠ 未分类 {len(unmapped)} 个特征 (需补 FAMILY_MAP):')
        for u in unmapped: print(f'    {u}')
        raise RuntimeError('请补全 FAMILY_MAP 后重跑')

    fmap = pd.DataFrame(rows)
    fmap = fmap.sort_values(['family_main', 'family_sub', 'short_label']).reset_index(drop=True)

    out_csv = PREP_DIR / 'metabolite_family_map.csv'
    atomic_write_csv(fmap, out_csv)
    print(f'\n✓ 映射写入: {out_csv.relative_to(ROOT)}')

    # 家族汇总
    print('\n=== 家族成员数 ===')
    for fm in fmap['family_main'].unique():
        sub = fmap[fmap['family_main'] == fm]
        n_in80 = int(sub['in_main_80'].sum())
        n_cand80 = int(sub['is_candidate_80'].sum())
        n_cand50 = int(sub['is_candidate_50'].sum())
        print(f'  {fm:30s} 总 {len(sub):3d}  在 80 轨 {n_in80:3d}  '
              f'候选 (80/50): {n_cand80}/{n_cand50}')

    # markdown summary
    md_lines = ['# 代谢物 → 家族映射汇总\n',
                '映射表: `data/02_preprocessed/metabolite_family_map.csv`  ',
                '生成: `scripts/12_metabolite_family_map.py`\n',
                '## 7 大家族总览\n',
                '| 家族 | 总数 | 在 80 主轨 | 80 候选 | 50 候选 |',
                '|---|---|---|---|---|']
    for fm in fmap['family_main'].unique():
        sub = fmap[fmap['family_main'] == fm]
        md_lines.append(f'| **{fm}** | {len(sub)} | {int(sub["in_main_80"].sum())} | '
                        f'{int(sub["is_candidate_80"].sum())} | '
                        f'{int(sub["is_candidate_50"].sum())} |')

    md_lines.append('\n## 各家族详细成员\n')
    for fm in fmap['family_main'].unique():
        sub = fmap[fmap['family_main'] == fm]
        md_lines.append(f'### {fm}  ({len(sub)} 个)\n')
        md_lines.append('| short | family_sub | substrate | enzyme | ev | non_enz | in_80 | cand_80 | cand_50 |')
        md_lines.append('|---|---|---|---|---|---|---|---|---|')
        for _, r in sub.iterrows():
            tag80 = '★' if r['is_candidate_80'] else ('✓' if r['in_main_80'] else '-')
            tag50 = '★' if r['is_candidate_50'] else '-'
            md_lines.append(f'| `{r["short_label"]}` | {r["family_sub"]} | '
                            f'{r["substrate"]} | {r["enzyme"]} | '
                            f'{r["enzyme_evidence_level"]} | {r["non_enzymatic_contribution"]} | '
                            f'{"✓" if r["in_main_80"] else "-"} | '
                            f'{tag80} | {tag50} |')
        md_lines.append('')

    md_lines.append('---\n**标记说明**:')
    md_lines.append('- ✓ = 在 80 主轨内; ★ = 差异代谢物候选 (p_limma<0.05 + |log2FC|≥log2(1.2) + 离群稳健)')
    md_lines.append('- **ev (enzyme_evidence_level)**:')
    md_lines.append('  - **5** = 结构唯一可达 (PG/TX/epoxide/diol/HHT 等; 自由基不可达终产物)')
    md_lines.append('  - **4** = 酶主导 ≥70%, 少量非酶 (12/15-HETE, 14-/17-/20-HDoHE, 16/18-HETE)')
    md_lines.append('  - **3** = 酶 vs 非酶 ~50:50 (oxoODE, 13-HpODE)')
    md_lines.append('  - **2** = 非酶主导 ~60-70% (**16-HDoHE 2026-05-19 重归此**, 4-HDoHE, 9-HODE, 11-HEPE, 8-HETrE)')
    md_lines.append('  - **1** = 几乎纯非酶 (11/8-HETE, 7/8/10/11/13-HDoHE)')
    md_lines.append('  - **NA** = 前体 PUFA, 不参与氧化分类')
    md_lines.append('- **non_enz**: none / minor / substantial / dominant / near-exclusive / NA')
    md_lines.append('')
    md_lines.append('## 16-HDoHE 重归属决策 (2026-05-19)')
    md_lines.append('原 enzyme="CYP-ω", family_sub="HDoHE-ω (DHA)" → 现 enzyme="Non-enzymatic (radical)", family_sub="HDoHE-non-enz (DHA)".')
    md_lines.append('依据:')
    md_lines.append('1. **文献**: VanRollins 2008 J Lipid Res — CYP4-DHA 主产物是 19/20/22-HDHA, **16-HDHA 不在 CYP4 偏好位置 (ω-7)**; Yin & Porter 2005 — DHA 自由基氧化 C15 双烯丙位 H 脱氢 → C16 radical addition → 16-OH + 17E 双键重排 (panel 命名 "16-Hydroxy-...,17E,19Z-" 完全契合该机制).')
    md_lines.append('2. **数据 (实证)**: 单分子 ANCOVA log2FC=+0.173 (p=0.19), 与非酶组 7/10/11/13-HDoHE 中位 (+0.17) 高度同质, 远低于纯 CYP 20-HDoHE (+0.333, p=0.022); 加入"(16+20)/DHA"比值反而稀释纯 CYP 信号 (FC 1.31× vs 单 20 的 1.35×); 加入"(7+10+11+13+16)/DHA"自氧化指数效应几乎不变.')

    md_path = TABLES / 'family_summary.md'
    md_path.write_text('\n'.join(md_lines), encoding='utf-8')
    print(f'✓ 摘要写入: {md_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
