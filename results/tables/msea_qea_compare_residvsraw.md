# MSEA QEA 本地结果对比: residualized vs raw

本地 GlobalTest QEA, n=165, K-pathways=5 (KEGG hsa 手工策划), 10000 perms.

## residualized (协变量 age + GA + year 已洗)

| pathway_id   | pathway_name                            |   pathway_size_K |    Q_obs |   p_perm |     q_bh |   frac_up |
|:-------------|:----------------------------------------|-----------------:|---------:|---------:|---------:|----------:|
| hsa00590     | Arachidonic acid metabolism             |               34 | 210.897  | 0.20158  | 0.587841 |     0.794 |
| hsa00592     | alpha-Linolenic acid metabolism         |               19 | 181.59   | 0.379362 | 0.587841 |     0.842 |
| hsa01040     | Biosynthesis of unsaturated fatty acids |                5 |  40.0403 | 0.455454 | 0.587841 |     0     |
| hsa04923     | Regulation of lipolysis in adipocytes   |                3 |  47.6692 | 0.530347 | 0.587841 |     0.333 |
| hsa00591     | Linoleic acid metabolism                |               13 | 100.102  | 0.587841 | 0.587841 |     0.385 |

## raw (未残差化 ComBat 后矩阵)

| pathway_id   | pathway_name                            |   pathway_size_K |    Q_obs |   p_perm |     q_bh |   frac_up |
|:-------------|:----------------------------------------|-----------------:|---------:|---------:|---------:|----------:|
| hsa00590     | Arachidonic acid metabolism             |               34 | 205.397  | 0.256674 | 0.675932 |     0.794 |
| hsa00591     | Linoleic acid metabolism                |               13 | 156.378  | 0.40016  | 0.675932 |     0.385 |
| hsa00592     | alpha-Linolenic acid metabolism         |               19 | 187.505  | 0.405559 | 0.675932 |     0.895 |
| hsa01040     | Biosynthesis of unsaturated fatty acids |                5 |  27.6312 | 0.60324  | 0.682232 |     0     |
| hsa04923     | Regulation of lipolysis in adipocytes   |                3 |  34.8304 | 0.682232 | 0.682232 |     0.333 |

## 解读
- residualized 是主结果 (与脚本 09 ANCOVA / 脚本 16 家族 GSEA 框架一致)
- raw 用作对照, 看协变量调整是否改变富集 p 值显著性顺序
- 与脚本 17 KEGG GSEA-rank (rank-based) 互校: 若两者方向一致, 信号更稳健