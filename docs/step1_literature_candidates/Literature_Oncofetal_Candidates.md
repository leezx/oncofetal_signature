---
title: Step 1 — Literature Oncofetal Candidates
date: 2026-10-01
status: draft
source_corpus: core_ref/markdown (55 files)
---

# Step 1 — Literature-curated candidate intestinal oncofetal markers

> [!info] 范围
> 本表只负责产生 **Literature Oncofetal Candidates**，**不定义 Core**。Core 的判定留给后续 Level 1/2（独立数据轴：fetal vs adult、CRC vs normal、cross-species、epithelial-intrinsic）。

## 方法

- 语料：`core_ref/markdown/` 共 55 个全文文件；逐篇全文阅读，提取文中**明确**作为 fetal / fetal-like-regenerative / oncofetal 状态 marker、签名成员、验证基因或功能检验基因的条目；只记录原文陈述，不做外部推断。
- 每条证据标注证据层级：
  - **A = true fetal intestine**：胎儿组织体内，或胎儿肠来源 organoid/spheroid
  - **B = fetal-like / regenerative**：成年非肿瘤组织（损伤、DSS、照射、蠕虫、YAP 激活、TGFB1 等）
  - **C = cancer / oncofetal**：CRC 或其他肿瘤
- 鼠源基因统一映射为人 HGNC 符号（Ctgf→CCN2，Cyr61→CCN1；Ly6a 无人类直系同源基因，单列）。
- 引用格式经 Crossref DOI 核对；综述文献标注 `(review)`，其引用的原始研究见 `evidence_table.tsv` 的 `primary_ref_cited_by_review` 列。
- 注释列格式：`[层级]` + 每个层级最具代表性的一条证据（优先原始研究、优先有功能实验）+ `功能证据`（如有）。注释为对原文的英文转述。

## 分层（为避免把通路/下游基因混进 marker 列表）

| Tier | 定义 | 基因数 |
|---|---|---|
| T1 | 肠道/CRC 语境中被直接作为 fetal / regenerative / oncofetal **状态 marker** 的基因 | 82 |
| T2 | 同一程序中的**调控因子/通路/过程基因**（YAP/TAZ/TEAD、AP-1、NF-κB、p53、SWI/SNF、机械感受、自噬、EMT-TF、Wnt 拮抗剂、生长因子等） | 73 |
| T3 | **非肠道** oncofetal/fetal-like 证据（HCC 内皮/巨噬细胞、胃癌基质、胃 SPEM、肝、肺、乳腺、胰腺），仅作参照 | 52 |
| Down | 在 fetal-like / oncofetal 状态中**下调**的基因（阴性 marker） | 51 |

各表内排序：证据层级覆盖数（A+B+C 优先）→ 是否有 A → 原始研究篇数。

> [!warning] 语料局限（影响 A 层级判断）
> - 全语料仅 **12 个基因**有上调方向的 A 层级证据（ANXA1, Ly6a, TACSTD2, GJA1, IL33, SPP1, SPRR1A, YAP1, SMARCA4, SMARCC1；以及肝脏的 PLVAP, FOLR2）。且多数来自 Yui 2018（胎儿 organoid）或综述（Clevers 2026, Viragova 2024），而非胎肠体内数据。
> - Task brief 中依赖的关键原始文献 **不在语料中**：Mustata 2013（fetal spheroid; TACSTD2/GJA1 体内 E14）、Fawkner-Corbett 2021 / Elmentaite 2021（human fetal gut atlas; CLU in fetal progenitors）、Pikkupeura 2023（fetal YAP: ANKRD1/CTGF/CYR61/CLU/TNFRSF12A）。因此 **CLU 在本表中只有 B/C 证据、TNFRSF12A 仅 1 篇 B 证据**——这是语料缺口，不是生物学结论。
> - 两个文件内容错误：`A patient-derived organoid model captures fetal-like plasticity in colorectal cancer.md`（实为 Meacham & Morrison 2013）与 `Presence of onco-fetal neighborhoods in hepatocellular carcinoma….md`（实为 Shani 2020 乳腺癌 IL-33）。两篇目标论文需重新下载。
> - Mzoughi 2025（51-gene OnF / huOnF）、Moorman 2025（113-gene fetal / 14-gene core）、Chen 2023（GSEA gene sets）的完整签名基因只在补充表中，本表只收录正文/图中点名的基因。完整签名属于 Level 0 benchmarking layer，建议单独抓取补充表。
> - 10 篇文献无符合标准的 marker（详见 `paper_inventory.tsv`）。

## 文件

- `evidence_table.tsv` — 逐条证据（563 行）：Gene | Paper | Species | Tissue | fetal_in_vivo | fetal_organoid | regeneration | CRC | other_cancer | functional_evidence | evidence_class | direction | provenance | note …
- `gene_summary.tsv` — 基因层面汇总（tier、A/B/C、各证据轴、原始研究/综述篇数）
- `paper_inventory.tsv` — 55 个文件的引用、DOI、贡献条目数、问题标注

## T1 — Intestinal oncofetal marker candidates

| Oncofetal marker | 参考文献 | 注释（为何入选） |
|---|---|---|
| ANXA1 | Yui et al., Cell Stem Cell, 2018; Ayyaz et al., Nature, 2019; Cañellas-Socias et al., Nature, 2022; Gil Vazquez et al., Cell Stem Cell, 2022; Álvarez-Varela et al., Nat Cancer, 2022; Chen et al., Cell Stem Cell, 2023; Qin et al., Cell, 2023; Hartl et al., Sci Adv, 2024; Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025; van der Net et al., Cell Rep, 2025; Buissant des Amorie et al., Nature, 2026; Sphyris et al., Cancers (Basel), 2021 (review); Tape, Trends Cancer, 2024 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: Fetal-enriched annexin used as a fetal-state marker (Yui 2018); B: Regenerative marker up in Trp53 cKO colon 2 mo post-DSS and in Wnt-independent Trp53 KO-WI organoids (YAP target) (Hartl 2024); C: Named fetal intestinal progenitor and YAP target marker up in chemo-resistant residual cells (Álvarez-Varela 2022) ‖ 功能证据: Trp53 cKO locks state; qPCR/bulk RNA (Hartl 2024); reporter design (Mzoughi 2025) |
| Ly6a (mouse; no human ortholog) | Nusse et al., Nature, 2018; Yui et al., Cell Stem Cell, 2018; Ayyaz et al., Nature, 2019; Cheung et al., Cell Stem Cell, 2020; Gil Vazquez et al., Cell Stem Cell, 2022; Chen et al., Cell Stem Cell, 2023; Hartl et al., Sci Adv, 2024; Mzoughi et al., Nat Genet, 2025; van der Net et al., Cell Rep, 2025; Hageman et al., Dev Cell, 2020 (review); Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review); Beumer & Clevers, Nat Rev Mol Cell Biol, 2021 (review); Sphyris et al., Cancers (Basel), 2021 (review); Palikuqi et al., Cold Spring Harb Perspect Biol, 2022 (review); Higa & Nakayama, Cancer Sci, 2024 (review); Tape, Trends Cancer, 2024 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: Sca1 marks primitive fetal epithelium in vivo and fetal organoids; multiple Ly6 members up (Yui 2018); B: Defining marker used to isolate the fetal-like repairing epithelium; YAP/TAZ-dependent (Yui 2018); C: Top oncofetal marker; locus used for OnF reporter construction (Mzoughi 2025) ‖ 功能证据: FACS sorting of EpCAM+Sca1high repairing epithelium; absent in YAP/TAZ cDKO cysts (Yui 2018); reporter design (Mzoughi 2025) ‖ 注意: Nusse 2018 报告 E15.5 胎肠上皮表达 Sca-1，与 task brief 引用的 Mustata 2013（E15 几乎不表达）观察相冲突；Mustata 原文不在本语料中，待核。无人类直系同源基因。 |
| TACSTD2 | Yui et al., Cell Stem Cell, 2018; Cañellas-Socias et al., Nature, 2022; Moorman et al., Nature, 2025; Ogden et al., Cell Genom, 2025; Sphyris et al., Cancers (Basel), 2021 (review); Goldenring & Mills, Gastroenterology, 2022 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: Fetal progenitor marker elevated in fetal organoids (Yui 2018); B: Fetal marker induced during repair and by mechano/YAP activation (Yui 2018); C: Injury-repair module member; protein higher in metastases than matched primaries (Moorman 2025) ‖ 功能证据: Rho-inhibition reduces Tacstd2 induction in ApcKO collagen organoids (Yui 2018) ‖ 本语料中 A 层级证据来自 Yui 2018（胎儿 organoid）及综述；Mustata 2013 原始 E14 体内证据不在语料中。 |
| GJA1 | Nusse et al., Nature, 2018; Mo et al., Hum Cell, 2023; Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: Gene described as typically enriched in organoids grown from fetal gut (Clevers 2026); B: Fetal marker up-regulated in Sca-1+ granuloma crypt cells in vivo (Nusse 2018); C: Defines oncofetal TEC subtype; poor prognosis (Mo 2023) |
| IL33 | Gregorieff et al., Nature, 2015; Burdziak et al., Science, 2023; Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: Gene described as typically enriched in organoids grown from fetal gut (Clevers 2026); B: Induced when Yap is made nuclear by Lats1/2 deletion (Gregorieff 2015); C: Epithelium-derived ligand of a feedback loop driving plastic pre-neoplastic states; cell_type: pancreatic epithelium (Burdziak 2023) ‖ 功能证据: shRNA knockdown in Kras-mutant epithelium (KC-shIl33) shifts epithelial and immune cell states and retains cells in progenitor states (Burdziak 2023) |
| SPP1 | Nusse et al., Nature, 2018; Burdziak et al., Science, 2023; Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: Gene described as typically enriched in organoids grown from fetal gut (Clevers 2026); B: Fetal marker up-regulated in Sca-1+ granuloma crypt cells in vivo (Nusse 2018); C: Co-marker of gastric (E6) plastic epithelial state; cell_type: pancreatic epithelium (Burdziak 2023) |
| SPRR1A | Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B] A: Gene described as typically enriched in organoids grown from fetal gut (Clevers 2026); B: Upregulated in crypts near worm-induced granulomas (fetal reversion), driven by interferon-gamma (Clevers 2026) ‖ 原文写作 Sprr1（未指明亚型），映射为 SPRR1A 存在不确定性。 |
| CCN2 | Yui et al., Cell Stem Cell, 2018; Cañellas-Socias et al., Nature, 2022; Heinz et al., Cancer Res, 2022; Álvarez-Varela et al., Nat Cancer, 2022; Chen et al., Cell Stem Cell, 2023; Hartl et al., Sci Adv, 2024; van der Net et al., Cell Rep, 2025; Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review) | [B·C] B: YAP/TAZ target gene marking the reprogrammed fetal-like state (Yui 2018); C: YAP target enriched in growth-stalled micro-organoids and earliest metastatic cells, the fetal-like YAP state (Heinz 2022) ‖ 功能证据: qPCR; verteporfin/Rho inhibition reduce induction; reversible on Matrigel transfer (Yui 2018) |
| CLU | Gregorieff et al., Nature, 2015; Ayyaz et al., Nature, 2019; Gil Vazquez et al., Cell Stem Cell, 2022; Chen et al., Cell Stem Cell, 2023; Qin et al., Cell, 2023; Hageman et al., Dev Cell, 2020 (review); Beumer & Clevers, Nat Rev Mol Cell Biol, 2021 (review); Sphyris et al., Cancers (Basel), 2021 (review); Higa & Nakayama, Cancer Sci, 2024 (review); Tape, Trends Cancer, 2024 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [B·C] B: Hallmark fetal/regenerative (revival) marker elevated in Ki67+ regenerating crypt cells post-IR and induced within 6h of TGFB1 in organoids (Chen 2023); C: revCSC state accessible in CRC-mutant epithelia; used as the revCSC axis of the stemness score (Qin 2023) ‖ 功能证据: Clu-CreERT2 lineage tracing yields crypt-villus ribbons and new LGR5+ CBCs; Clu-DTA ablation impairs regeneration and DSS survival (Ayyaz 2019); TGFB1 sufficient to induce in organoids; induction lost in Tgfbr2KO and with TGFBR inhibitors; scRNA/RNA velocity (Chen 2023) |
| L1CAM | Damelin et al., Cancer Res, 2011; Ganesh et al., Nat Cancer, 2020; Cañellas-Socias et al., Nature, 2022; Moorman et al., Nature, 2025; De Angelis et al., Cancers (Basel), 2021 (review); Sphyris et al., Cancers (Basel), 2021 (review); Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Tape, Trends Cancer, 2024 (review) | [B·C] B: Induced during post-colitis regeneration and in normal colon organoid formation; required for epithelial repair (Ganesh 2020); C: Marks metastasis-initiating, chemoresistant organoid-forming cells; absent in normal colon, enriched at invasion front and metastases (Ganesh 2020) ‖ 功能证据: Epithelial and Lgr5-lineage L1cam deletion impairs healing and survival after DSS colitis; IHC in regenerating crypts (Ganesh 2020); CRISPR KO and inducible shRNA abolish organoid regeneration, orthotopic growth and liver metastasis; IHC/flow (Ganesh 2020) |
| CD44 | Damelin et al., Cancer Res, 2011; Ganesh et al., Nat Cancer, 2020; Chen et al., Cell Stem Cell, 2023; Hartl et al., Sci Adv, 2024; Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review); De Angelis et al., Cancers (Basel), 2021 (review); Goldenring & Mills, Gastroenterology, 2022 (review); Viragova et al., Cell Stem Cell, 2024 (review) | [B·C] B: Stem/proliferative marker in regenerating clusters post-IR; regeneration marker with TGFB1-induced chromatin accessibility (Chen 2023); C: Higher on L1CAM-high than L1CAM-low tumour cells (Ganesh 2020) |
| CCN1 | Gregorieff et al., Nature, 2015; Yui et al., Cell Stem Cell, 2018; Cañellas-Socias et al., Nature, 2022; Heinz et al., Cancer Res, 2022; Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review) | [B·C] B: YAP/TAZ target gene induced with the fetal-like reprogramming (Yui 2018); C: YAP target gene of the micro-organoid/micrometastatic fetal-like state (Heinz 2022) ‖ 功能证据: qPCR; mechano-inhibition lowers expression (Yui 2018) |
| ANKRD1 | Yui et al., Cell Stem Cell, 2018; Heinz et al., Cancer Res, 2022; Hartl et al., Sci Adv, 2024; van der Net et al., Cell Rep, 2025 | [B·C] B: YAP target gene massively up in locked-regenerative Trp53 KO-WI organoids (Hartl 2024); C: YAP target enriched in micro-organoids and cluster 6 (day 1) metastatic cells (Heinz 2022) |
| SOX9 | Messier et al., J Cell Physiol, 2016; Chen et al., Cell Stem Cell, 2023; Qin et al., Cell, 2023; Goldenring & Mills, Gastroenterology, 2022 (review); Viragova et al., Cell Stem Cell, 2024 (review) | [B·C] B: Pro-regenerative TF induced by IR and TGFB1, co-expressed with Tgfbr2/Clu; required for fetal/regenerative response (Chen 2023); C: Example tumor suppressor in APC pathway bivalently marked (H3K4me3+H3K27me3) in hESC, MCF7 and patient tumors; 'oncofetal epigenetic' signature; direction=bivalent, not expression (Messier 2016) ‖ 功能证据: Sox9KO blocks TGFB1-induced spheroids/fetal genes and reduces OLFM4+ regenerative crypts after IR; SMAD4 ChIP binding (Chen 2023) |
| ITGB1 | Chen et al., Cell Stem Cell, 2023; van der Net et al., Cell Rep, 2025; Centonze et al., Cancer Discov, 2026 | [B·C] B: Adhesion receptor upregulated during IR regeneration and by TGFB1; part of YAP-overlapping wound-healing program (Chen 2023); C: Collagen I-binding integrin mediating fetal-like transition (van der Net 2025) |
| MEX3A | Cañellas-Socias et al., Nature, 2022; Álvarez-Varela et al., Nat Cancer, 2022; Morgan et al., Br J Cancer, 2018 (review); Cao et al., MedComm, 2023 (review); Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Higa & Nakayama, Cancer Sci, 2024 (review); Tape, Trends Cancer, 2024 (review) | [B·C] B: Review cites Mex3a as marker of slow-cycling Lgr5+ cells that survive injury and chemo/radiation (Morgan 2018); C: Marks latent slow-cycling LGR5+ persister cells that convert to a YAP+ revival/fetal-like state after chemotherapy (Álvarez-Varela 2022) ‖ 功能证据: lineage tracing with Mex3a-CreERT2 regenerates organoids and metastases after FOLFIRI; Mex3a KO organoids chemosensitive; iCaspase9 ablation plus FOLFIRI prolongs response (Álvarez-Varela 2022) |
| BASP1 | Ayyaz et al., Nature, 2019; Álvarez-Varela et al., Nat Cancer, 2022; Sphyris et al., Cancers (Basel), 2021 (review) | [B·C] B: Listed in the quiescent SSC2c revival stem cell signature (Ayyaz 2019); C: Fetal progenitor / revival stem cell marker up during FOLFIRI and down after recovery (Álvarez-Varela 2022) |
| MSLN | Gregorieff et al., Nature, 2015; Centonze et al., Cancer Discov, 2026; Viragova et al., Cell Stem Cell, 2024 (review) | [B·C] B: Induced when Yap is made nuclear by Lats1/2 deletion (Gregorieff 2015); C: Listed as marker of KRAS-activated HRCs overlapping basal PDAC program (Centonze 2026) |
| AMOTL2 | Cheung et al., Cell Stem Cell, 2020; Heinz et al., Cancer Res, 2022 | [B·C] B: YAP target marking wound-healing state (Cheung 2020); C: Known YAP target enriched in the earliest (day 1) metastatic cell cluster (Heinz 2022) |
| CLDN4 | Cañellas-Socias et al., Nature, 2022; Qin et al., Cell, 2023 | [B·C] B: Epithelial progenitor gene high in fibroblast-induced revCSC cluster (Qin 2023); C: Tight junction component upregulated in the core HRC program (Cañellas-Socias 2022) |
| KLF6 | Cheung et al., Cell Stem Cell, 2020 | [B·C] B: Marks YAP-induced Wnt-low wound-healing cell state absent in homeostasis (Cheung 2020); C: Klf6 up as YAP reprograms Apc-mutant stem cells (Cheung 2020) |
| KRT19 | Morgan et al., Br J Cancer, 2018 (review); De Angelis et al., Cancers (Basel), 2021 (review) | [B·C] B: Review cites KRT19+/Lgr5-negative cells as responsible for cancer tissue regrowth after irradiation (De Angelis 2021); C: KRT19+ Lgr5-negative cells regrow cancer tissue after irradiation; stress-induced stem function (De Angelis 2021) |
| EMP1 | Cañellas-Socias et al., Nature, 2022; Moorman et al., Nature, 2025; Ogden et al., Cell Genom, 2025; van der Net et al., Cell Rep, 2025; Buissant des Amorie et al., Nature, 2026; Centonze et al., Cancer Discov, 2026; Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Tape, Trends Cancer, 2024 (review) | [C] C: Marker gene chosen for High Relapse Cells; marks invasion fronts, tumor buds and micrometastases (Cañellas-Socias 2022) ‖ 功能证据: CRISPR mNeonGreen knock-in reporter at EMP1 locus used to screen oncofetal-state inducers (Buissant des Amorie 2026); Emp1-iCaspase9 ablation prevents metastatic relapse; knock-in reporter tracking; RNA-FISH in patient CRC (Cañellas-Socias 2022) |
| LAMC2 | Cañellas-Socias et al., Nature, 2022; Ogden et al., Cell Genom, 2025; Buissant des Amorie et al., Nature, 2026; Centonze et al., Cancer Discov, 2026 | [C] C: coreHRC and partial-EMT gene; known tumor budding marker correlating with EMP1 in patient samples (Cañellas-Socias 2022) ‖ 功能证据: YAP knockdown/TEADi did not change Lamc2, showing YAP-independence (Cañellas-Socias 2022) |
| PLAUR | Cañellas-Socias et al., Nature, 2022; Gil Vazquez et al., Cell Stem Cell, 2022; Ogden et al., Cell Genom, 2025; Centonze et al., Cancer Discov, 2026 | [C] C: coreHRC partial-EMT gene; protease receptor associated with invasion and extravasation (Cañellas-Socias 2022) |
| ITGA2 | Cañellas-Socias et al., Nature, 2022; van der Net et al., Cell Rep, 2025; Centonze et al., Cancer Discov, 2026 | [C] C: Collagen I-binding integrin mediating fetal-like transition (van der Net 2025) ‖ 功能证据: BTT-3033 blocking reduces Sca-1high cells; with TRPV4 inhibition fully blocks (van der Net 2025) |
| PROX1 | Moorman et al., Nature, 2025; Ogden et al., Cell Genom, 2025; Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Higa & Nakayama, Cancer Sci, 2024 (review); Tape, Trends Cancer, 2024 (review) | [C] C: Fetal-state TF repressing non-intestinal lineages; marks CDX2low invasion-front cells (Moorman 2025) ‖ 功能证据: shRNA knockdown immunofluorescence (Moorman 2025) |
| ANXA2 | Cheung et al., Cell Stem Cell, 2020; Centonze et al., Cancer Discov, 2026 | [C] C: HRC/fetal-like marker gene co-expressed with Emp1 in metastatic clusters (Centonze 2026) |
| ITGB4 | Cañellas-Socias et al., Nature, 2022; Centonze et al., Cancer Discov, 2026 | [C] C: Integrin enriched in HRC differential expression analysis (Cañellas-Socias 2022) |
| JUP | Cañellas-Socias et al., Nature, 2022; Centonze et al., Cancer Discov, 2026 | [C] C: Junction component of core HRC program; linked to clustered tumor cell dissemination (Cañellas-Socias 2022) |
| LAMA3 | Cañellas-Socias et al., Nature, 2022; Centonze et al., Cancer Discov, 2026 | [C] C: coreHRC gene assigned to the partial EMT module upregulated in HRCs (Cañellas-Socias 2022) |
| S100A6 | Cheung et al., Cell Stem Cell, 2020; Mzoughi et al., Nat Genet, 2025 | [C] C: Listed as a top oncofetal marker; its locus used to build OnF sLCR reporter (Mzoughi 2025) ‖ 功能证据: reporter design (Mzoughi 2025) |
| VEGFA | Sharma et al., Cell, 2020; Ogden et al., Cell Genom, 2025 | [C] C: Tumor-hepatocyte ligand predicted to induce PLVAP in endothelium; epithelial/tumor-cell derived (Sharma 2020) |
| IL1RN | Gregorieff et al., Nature, 2015; Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [B] B: Yap-induced regenerative/inflammatory gene (Gregorieff 2015) |
| CXADR | Ayyaz et al., Nature, 2019; Sphyris et al., Cancers (Basel), 2021 (review) | [B] B: Listed in the quiescent SSC2c revival stem cell signature (Ayyaz 2019) |
| PROM1 | Ganesh et al., Nat Cancer, 2020; De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Higher on L1CAM-high than L1CAM-low tumour cells (Ganesh 2020) |
| ANXA3 | Álvarez-Varela et al., Nat Cancer, 2022 | [C] C: Named fetal progenitor / revival stem cell marker tracked during FOLFIRI and recovery (Álvarez-Varela 2022) |
| CD274 | Cañellas-Socias et al., Nature, 2022 | [C] C: Immunomodulatory gene elevated in undifferentiated HRCs of micrometastases (Cañellas-Socias 2022) |
| CD70 | Moorman et al., Nature, 2025 | [C] C: Injury-repair module gene linked to regeneration and therapy resistance (Moorman 2025) |
| COL1A1 | Yui et al., Cell Stem Cell, 2018 | [B] B: ECM gene strongly elevated in the repairing (fetal-like) epithelium itself (Yui 2018) |
| CXCL2 | Ogden et al., Cell Genom, 2025 | [C] C: Chemokine marking the hybrid intermediate (i)REC/stem state and a predicted iREC-inducing ligand (Ogden 2025) |
| CXCL3 | Ogden et al., Cell Genom, 2025 | [C] C: Chemokine upregulated in the intermediate REC/stem hybrid cancer state (Ogden 2025) |
| CYP2W1 | Moorman et al., Nature, 2025 | [C] C: Embryonic developmental gene elevated in treatment-naive tumour cells above normal stem cells (Moorman 2025) |
| DDR1 | Cañellas-Socias et al., Nature, 2022 | [C] C: Collagen-sensing receptor enriched in HRCs vs non-HRCs (Cañellas-Socias 2022) |
| DSC2 | Cañellas-Socias et al., Nature, 2022 | [C] C: Desmosome component upregulated in the core HRC program (Cañellas-Socias 2022) |
| EPHB2 | Ganesh et al., Nat Cancer, 2020 | [C] C: Higher on L1CAM-high than L1CAM-low tumour cells (Ganesh 2020) |
| EZR | Centonze et al., Cancer Discov, 2026 | [C] C: Cytoskeleton organizer listed among EpiHR HRC genes (Centonze 2026) |
| F11R | Centonze et al., Cancer Discov, 2026 | [C] C: Tight-junction gene listed in EpiHR HRC signature (Centonze 2026) |
| GPC1 | Moorman et al., Nature, 2025 | [C] C: Fetal-state readout rising when intestinal growth factors removed from metastasis organoids (Moorman 2025) |
| IDO1 | Cañellas-Socias et al., Nature, 2022 | [C] C: Immune-suppressive gene elevated in micrometastatic HRCs (Cañellas-Socias 2022) |
| IGFBP3 | Hartl et al., Sci Adv, 2024 | [B] B: Regenerative marker increased in locked-regenerative Trp53 cKO colon (Hartl 2024) |
| ITGA3 | Cañellas-Socias et al., Nature, 2022 | [C] C: Integrin enriched in HRC differential expression analysis (Cañellas-Socias 2022) |
| KLK8 | Centonze et al., Cancer Discov, 2026 | [C] C: Keratinization mediator in Emp1+ HRC subset (Centonze 2026) |
| KRT17 | Cañellas-Socias et al., Nature, 2022 | [C] C: Basal pancreatic-subtype marker that labels EMP1-high invasion fronts and tumor buds (Cañellas-Socias 2022) ‖ 功能证据: immunostaining; gained in MTO co-cultured with colon fibroblasts (Cañellas-Socias 2022) |
| KRT27 | Centonze et al., Cancer Discov, 2026 | [C] C: Keratinization gene in Emp1+ HRC subset, KRAS-driven (Centonze 2026) |
| KRT5 | Centonze et al., Cancer Discov, 2026 | [C] C: Basal keratin upregulated in MAPK/KRAS-mutant CRC alongside HRC program (Centonze 2026) |
| KRT6A | Centonze et al., Cancer Discov, 2026 | [C] C: Basal keratin upregulated in MAPK/KRAS-mutant CRC alongside HRC program (Centonze 2026) |
| KRT79 | Centonze et al., Cancer Discov, 2026 | [C] C: Keratinization gene in Emp1+ HRC subset, KRAS-driven (Centonze 2026) |
| KRT80 | Centonze et al., Cancer Discov, 2026 | [C] C: Keratinization gene in Emp1+ HRC subset, KRAS-driven (Centonze 2026) |
| LAMB3 | Centonze et al., Cancer Discov, 2026 | [C] C: Hemidesmosome component of EpiHR high-relapse-cell signature (Centonze 2026) |
| Ly6c1 (mouse) | Gregorieff et al., Nature, 2015 | [B] B: Yap-dependent regenerative gene; Ly6c1 has no direct human ortholog (Gregorieff 2015) |
| MMP7 | Buissant des Amorie et al., Nature, 2026 | [C] C: Top invasive-front gene, named as regenerative/fetal state marker (Buissant des Amorie 2026) |
| OSMR | Moorman et al., Nature, 2025 | [C] C: Injury-repair module gene linked to regeneration and therapy resistance (Moorman 2025) |
| PCDH1 | Cañellas-Socias et al., Nature, 2022 | [C] C: Adhesion/junction component upregulated in the core HRC program (Cañellas-Socias 2022) |
| PERP | Centonze et al., Cancer Discov, 2026 | [C] C: Keratinization mediator in Emp1+ HRC subset (Centonze 2026) |
| PKM | Hartl et al., Sci Adv, 2024 | [B] B: Rate-limiting glycolytic enzyme up in locked regenerative state; noted as typical of fetal tissues and cancers (Hartl 2024) ‖ 功能证据: proteomics, IF; 2-DG glycolysis inhibition blocks KO-WI growth (Hartl 2024) |
| PTK7 | Moorman et al., Nature, 2025 | [C] C: Member of 14-gene core fetal progenitor set shared by four patients (Moorman 2025) |
| Reg3b (mouse) | Yui et al., Cell Stem Cell, 2018 | [B] B: C-type lectin co-induced with Sca1 in the repairing epithelium (Yui 2018) |
| REG3G | Yui et al., Cell Stem Cell, 2018 | [B] B: Co-induced with Sca1 in the repairing epithelium (Yui 2018) |
| RHOF | Centonze et al., Cancer Discov, 2026 | [C] C: Cytoskeleton organizer listed among EpiHR HRC genes (Centonze 2026) |
| RRAS | Centonze et al., Cancer Discov, 2026 | [C] C: Actin cytoskeleton organizer listed in EpiHR HRC signature (Centonze 2026) |
| S100A11 | Cheung et al., Cell Stem Cell, 2020 | [C] C: Injury signature gene upregulated during YAP reprogramming (Cheung 2020) |
| TGM1 | Centonze et al., Cancer Discov, 2026 | [C] C: Keratinization mediator in Emp1+ HRC subset (Centonze 2026) |
| TJP3 | Centonze et al., Cancer Discov, 2026 | [C] C: Tight-junction gene listed in EpiHR HRC signature (Centonze 2026) |
| TMEM132A | Moorman et al., Nature, 2025 | [C] C: Fetal-state readout rising when intestinal growth factors removed from metastasis organoids (Moorman 2025) |
| TNFRSF12A | Gregorieff et al., Nature, 2015 | [B] B: Yap-induced regenerative gene co-expressed with Olfm4 in regenerating stem cells (Gregorieff 2015) ‖ 注意: 本语料中仅 Gregorieff 2015 一篇（B 层级）；task brief 提到的 Pikkupeura fetal YAP 证据不在语料中。 |
| KRT20 | Morgan et al., Br J Cancer, 2018 (review); De Angelis et al., Cancers (Basel), 2021 (review); Higa & Nakayama, Cancer Sci, 2024 (review) | [C] C: Lgr5-negative KRT20+ cells drive tumor regrowth and regenerate Lgr5+ cells after ablation (De Angelis 2021) |
| ALCAM | De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Review: used with CD44 for cCSC isolation; downregulated in infiltrating edge cells (De Angelis 2021) |
| AXL | Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review) | [C] C: YAP/TAZ target correlated with poor survival and part of the CMS4 gene set (Mouillet-Richard & Laurent-Puig 2020) |
| HOPX | Walter et al., Cancers (Basel), 2021 (review) | [B] B: Named as marker of colitis-associated regenerative stem cells; review notes possible relevance to CRC (Walter 2021) |
| TDGF1 | De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Review (authors' own work): fluctuating marker reflecting plastic cCSC state (De Angelis 2021) |
| TFF3 | Goldenring & Mills, Gastroenterology, 2022 (review) | [C] C: Intestinal phenotype marker emerging after SPEM in Kras-driven metaplasia preceding dysplasia (Goldenring & Mills 2022) |

## T2 — Program regulators / pathway & process genes

| Oncofetal marker | 参考文献 | 注释（为何入选） |
|---|---|---|
| YAP1 | Gregorieff et al., Nature, 2015; Yui et al., Cell Stem Cell, 2018; Ayyaz et al., Nature, 2019; Cheung et al., Cell Stem Cell, 2020; Cañellas-Socias et al., Nature, 2022; Chen et al., Cell Stem Cell, 2023; Hartl et al., Sci Adv, 2024; Mzoughi et al., Nat Genet, 2025; van der Net et al., Cell Rep, 2025; Hageman et al., Dev Cell, 2020 (review); Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review); Beumer & Clevers, Nat Rev Mol Cell Biol, 2021 (review); Sphyris et al., Cancers (Basel), 2021 (review); Goldenring & Mills, Gastroenterology, 2022 (review); Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Higa & Nakayama, Cancer Sci, 2024 (review); Tape, Trends Cancer, 2024 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] A: YAP activation described as required to maintain fetal epithelial state and transiently in intestinal development (Clevers 2026); B: Active YAP marks regenerative state after DSS and persists in Trp53 cKO/KO-WI; YAP UP signature enriched (Hartl 2024); C: Nuclear Yap drives the regenerative program in Apc-mutant tumour-initiating cells (Gregorieff 2015) ‖ 功能证据: Yap1 epithelial knockout abolishes Clu+ cells and impairs regeneration; Lats1/2 knockout induces ectopic Clu-GFP+ cells (Ayyaz 2019); YAP S127A/5SA overexpression reprograms ISCs; lineage tracing shows loss of stemness (Cheung 2020) |
| WWTR1 | Yui et al., Cell Stem Cell, 2018; Cheung et al., Cell Stem Cell, 2020; Hageman et al., Dev Cell, 2020 (review); Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review); Sphyris et al., Cancers (Basel), 2021 (review); Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Higa & Nakayama, Cancer Sci, 2024 (review) | [B·C] B: TAZ required with YAP for establishing the fetal-like repairing state (Yui 2018); C: Deletion with Yap promotes colon tumor growth (Cheung 2020) ‖ 功能证据: YAP/TAZ cDKO blocks collagen growth, repair and Sca1 induction (Yui 2018); Yap/Taz deletion increases tumor burden (Cheung 2020) |
| PTK2 | van der Net et al., Cell Rep, 2025; Hageman et al., Dev Cell, 2020 (review); Sphyris et al., Cancers (Basel), 2021 (review) | [B·C] B: Rises with collagen 1 and integrin after injury; links ECM sensing to YAP-dependent regeneration (Hageman 2020); C: Active FAK elevated in collagen I; required for integrin-driven reprogramming (van der Net 2025) ‖ 功能证据: PF-573228 FAK inhibition reduces Sca-1high cells (van der Net 2025) |
| TP53 | Hartl et al., Sci Adv, 2024; Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review) | [B·C] B: p53 induced in regenerative epithelium after DSS; required to terminate fetal-like regenerative state (Hartl 2024); C: Chemotherapy-induced quiescence tied to the fetal-like state arises in TP53-wild-type tumour cells (Pérez-González 2023) ‖ 功能证据: cKO prevents exit from regenerative state; Nutlin-3a activation suppresses Wnt (Hartl 2024) |
| ATF3 | Mzoughi et al., Nat Genet, 2025; Goldenring & Mills, Gastroenterology, 2022 (review) | [B·C] B: Early-response transcription factor initiating the regenerative reprogramming program (Goldenring & Mills 2022); C: AP-1 subunit induced by active YAP in wild-type organoids (Mzoughi 2025) |
| NOTUM | Ogden et al., Cell Genom, 2025; Ramadan et al., Trends Cancer, 2022 (review) | [B·C] B: Review: Wnt antagonist secreted by mutant and aged cells that differentiates neighbouring normal ISCs (Ramadan 2022); C: Defines a cancer-specific stem NOTUM state (WNT antagonist), distinct from but related to non-canonical states (Ogden 2025) |
| SMARCA4 | Viragova et al., Cell Stem Cell, 2024 (review) | [A] A: SWI/SNF member required to maintain the fetal gene program (Viragova 2024) |
| SMARCC1 | Viragova et al., Cell Stem Cell, 2024 (review) | [A] A: SWI/SNF member gatekeeping fetal versus mature identity (Viragova 2024) |
| AREG | Gregorieff et al., Nature, 2015; Chen et al., Cell Stem Cell, 2023; Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [B] B: YAP target upregulated in epithelium during IR regeneration and upon TGFB1 treatment of organoids (Chen 2023) |
| EREG | Gregorieff et al., Nature, 2015; Yui et al., Cell Stem Cell, 2018; Beumer & Clevers, Nat Rev Mol Cell Biol, 2021 (review); Sphyris et al., Cancers (Basel), 2021 (review) | [B] B: Key Egfr ligand of the Yap regenerative program; functionally sufficient to rescue Yap loss (Gregorieff 2015) |
| EDN1 | Gregorieff et al., Nature, 2015; Chen et al., Cell Stem Cell, 2023 | [B] B: Regenerative gene induced in stem cells post-irradiation in a Yap-dependent manner (Gregorieff 2015) |
| FOS | Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025 | [C] C: AP-1 subunit reinforcing oncofetal state and lineage plasticity (Mzoughi 2025) ‖ 功能证据: knockdown dominant-negative overexpression (Mzoughi 2025) |
| FOSL1 | Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025 | [C] C: AP-1 TF with motif accessibility correlated to expression in (i)RECs (Ogden 2025) |
| FOSL2 | Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025 | [C] C: AP-1 subunit induced by active YAP in wild-type organoids (Mzoughi 2025) |
| JUND | Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025 | [C] C: AP-1 family footprint enriched in tumoroid-specific accessible regions (Mzoughi 2025) |
| TEAD1 | Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025 | [C] C: TEAD motif activity broad across transitional states of stem-oncofetal spectrum (Mzoughi 2025) |
| ZEB1 | Burdziak et al., Science, 2023; Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review); De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Marker of the EMT-like epithelial state captured in the tumorigenesis atlas; cell_type: pancreatic epithelium (Burdziak 2023) |
| TCF7 | Moorman et al., Nature, 2025; Higa & Nakayama, Cancer Sci, 2024 (review) | [C] C: Member of 14-gene core fetal progenitor set shared by four patients (Moorman 2025) |
| APCDD1 | Ogden et al., Cell Genom, 2025 | [C] C: WNT antagonist marking the cancer-specific stem NOTUM state (Ogden 2025) |
| ATG13 | Rehman et al., Cell, 2021 | [C] C: Autophagy component upregulated in diapause-like persister cells (Rehman 2021) |
| ATG16L2 | Rehman et al., Cell, 2021 | [C] C: Autophagy gene induced in irinotecan-induced persister cultures (Rehman 2021) |
| ATG2A | Rehman et al., Cell, 2021 | [C] C: Autophagy gene induced in irinotecan-induced persister cultures (Rehman 2021) |
| AXIN2 | Hartl et al., Sci Adv, 2024 | [B] B: Wnt target up after injury; remains elevated in Trp53 cKO and KO-WI organoids (Hartl 2024) |
| BATF2 | Mzoughi et al., Nat Genet, 2025 | [C] C: AP-1 subunit induced by active YAP in wild-type organoids (Mzoughi 2025) |
| BMP4 | Moorman et al., Nature, 2025 | [C] C: Endoderm development module gene; module most correlated with core fetal signature (Moorman 2025) |
| BMP7 | Moorman et al., Nature, 2025 | [C] C: Embryonic developmental gene elevated in treatment-naive tumour cells above normal stem cells (Moorman 2025) |
| CCND1 | Ayyaz et al., Nature, 2019 | [B] B: Gene distinguishing the irradiation-specific cluster 18 (SSC2) (Ayyaz 2019) |
| CDKN1A | Hartl et al., Sci Adv, 2024 | [B] B: p53 target induced in nearly all nuclei of regenerative crypts and Wnt-boosted organoids (Hartl 2024) |
| FOSB | Ogden et al., Cell Genom, 2025 | [C] C: AP-1 member whose regulon defines the (i)REC state; EMP1 is in its regulon (Ogden 2025) |
| HBEGF | Gregorieff et al., Nature, 2015 | [B] B: Egfr ligand listed in the Yap regenerative module (Gregorieff 2015) ‖ 功能证据: RNA-seq of Yap-deficient vs Yap-overexpressing organoids (Gregorieff 2015) |
| JUN | Mzoughi et al., Nat Genet, 2025 | [C] C: AP-1 family footprint activity rising with tumor progression (Mzoughi 2025) |
| JUNB | Ogden et al., Cell Genom, 2025 | [C] C: AP-1 member regulating (i)REC state genes (Ogden 2025) |
| MAFK | Mzoughi et al., Nat Genet, 2025 | [C] C: Motif activity increasing along malignancy continuum with AP-1 partners (Mzoughi 2025) |
| MAP1LC3B | Rehman et al., Cell, 2021 | [C] C: Autophagy marker increased in diapause-like persister cells (Rehman 2021) |
| MIF | Ayyaz et al., Nature, 2019 | [B] B: Gene distinguishing the irradiation-specific cluster 18 (SSC2) (Ayyaz 2019) |
| MKI67 | Hartl et al., Sci Adv, 2024 | [B] B: Proliferation increased in regenerative/locked-regenerative epithelium (Hartl 2024) |
| NFKB1 | Ogden et al., Cell Genom, 2025 | [C] C: NF-kB subunit identified as co-regulator of the inflammatory regenerative state (Ogden 2025) |
| NFKB2 | Ogden et al., Cell Genom, 2025 | [C] C: NF-kB subunit enriched in (i)REC states (Ogden 2025) |
| NKD1 | Ogden et al., Cell Genom, 2025 | [C] C: WNT antagonist marking the cancer-specific stem NOTUM state (Ogden 2025) |
| NPM1 | Ayyaz et al., Nature, 2019 | [B] B: Gene distinguishing the irradiation-specific cluster 18 (SSC2) (Ayyaz 2019) |
| OPTN | Rehman et al., Cell, 2021 | [C] C: Autophagy gene induced in irinotecan-induced persister cultures (Rehman 2021) |
| PIEZO1 | van der Net et al., Cell Rep, 2025 | [C] C: Ectopic activation induces reprogramming; not required under collagen I (van der Net 2025) |
| PTGS2 | Chen et al., Cell Stem Cell, 2023 | [B] B: YAP-related gene induced by TGFB1 in organoids; also pro-regenerative mesenchymal transcript (Chen 2023) ‖ 功能证据: TGFB1 induction blocked in Tgfbr2KO / TGFBR inhibitor (Chen 2023) |
| RELB | Ogden et al., Cell Genom, 2025 | [C] C: NF-kB subunit regulon defining the iREC state (Ogden 2025) |
| TEAD3 | Mzoughi et al., Nat Genet, 2025 | [C] C: TEAD footprint hyperactivated by RXR inhibition driving oncofetal conversion (Mzoughi 2025) |
| TEAD4 | Ogden et al., Cell Genom, 2025 | [C] C: TEAD TF enriched in the intermediate REC/stem state (Ogden 2025) |
| TGFBR2 | Chen et al., Cell Stem Cell, 2023 | [B] B: Receptor elevated day 3 post-IR specifically in Clu+/Ly6a+ regenerative epithelial cells; required for regeneration (Chen 2023) ‖ 功能证据: Epithelial Tgfbr2 KO restricts regeneration and abolishes TGFB1-induced fetal genes (Chen 2023) |
| TRPV4 | van der Net et al., Cell Rep, 2025 | [C] C: Mechanosensitive channel transducing collagen I forces to fetal-like reprogramming (van der Net 2025) ‖ 功能证据: RN-1734/GSK205 inhibition reduces Sca-1high emergence and Ca2+ influx (van der Net 2025) |
| ULK1 | Rehman et al., Cell, 2021 | [C] C: Autophagy driver of diapause-like drug-tolerant persister state; functionally required (Rehman 2021) ‖ 功能证据: ULK1 inhibitor plus irinotecan blocks persister state and kills cells (Rehman 2021) |
| VIM | Damelin et al., Cancer Res, 2011 | [C] C: EMT marker co-expressed with 5T4 in undifferentiated state (Damelin 2011) |
| WDR45 | Rehman et al., Cell, 2021 | [C] C: Autophagy gene induced in irinotecan-induced persister cultures (Rehman 2021) |
| WIPI1 | Rehman et al., Cell, 2021 | [C] C: Autophagy gene induced in irinotecan-induced persister cultures (Rehman 2021) |
| WIPI2 | Rehman et al., Cell, 2021 | [C] C: Autophagy gene induced in irinotecan-induced persister cultures (Rehman 2021) |
| WNT5A | Chen et al., Cell Stem Cell, 2023 | [B] B: YAP-related regeneration gene induced by TGFB1 in organoid epithelium; also induced in TGFB1-treated mesenchyme (Chen 2023) ‖ 功能证据: TGFB1 induction blocked in Tgfbr2KO / TGFBR inhibitor; ATAC open chromatin (Chen 2023) |
| WNT5B | Moorman et al., Nature, 2025 | [C] C: Endoderm development module gene; module most correlated with core fetal signature (Moorman 2025) |
| CDKN1B | Pérez-González et al., Nat Cancer, 2023 (review); Higa & Nakayama, Cancer Sci, 2024 (review) | [C] C: Marks dormant CSCs whose awakening depends on FAK-YAP activation (Pérez-González 2023) |
| ABCB1 | Cao et al., MedComm, 2023 (review) | [C] C: Framed as embryo/fetal-barrier transporter reactivated in drug-resistant cancer cells; not intestine-specific (Cao 2023) |
| ABCC1 | Cao et al., MedComm, 2023 (review) | [C] C: Framed as embryo/fetal-barrier transporter reactivated in drug-resistant cancer cells; not intestine-specific (Cao 2023) |
| ABCG2 | Cao et al., MedComm, 2023 (review) | [C] C: Framed as embryo/fetal-barrier transporter reactivated in drug-resistant cancer cells; not intestine-specific (Cao 2023) |
| ASCL2 | Ramadan et al., Trends Cancer, 2022 (review) | [B] B: Borderline: homeostatic ISC TF required for regenerative dedifferentiation; no fetal framing (Ramadan 2022) |
| CDH2 | Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review) | [C] C: YAP/TAZ-controlled EMT gene specifying the CMS4 phenotype (Mouillet-Richard & Laurent-Puig 2020) |
| CDKN1C | Higa & Nakayama, Cancer Sci, 2024 (review) | [C] C: Marks quiescent Lgr5+ CSCs activated by chemotherapy; also marks injury-induced intestinal stem cells (Higa & Nakayama 2024) |
| DKK2 | Ramadan et al., Trends Cancer, 2022 (review) | [B] B: Review: Wnt antagonist secreted by Apc-mutant cells against neighbouring normal ISCs (Ramadan 2022) |
| POLR1A | Ramadan et al., Trends Cancer, 2022 (review) | [C] C: Review: biosynthetic marker of functional CSCs in both LGR5+ and LGR5-negative cells near stroma (Ramadan 2022) |
| PRNP | Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review) | [C] C: Enriched in CMS4 tumours and required to sustain the YAP/TAZ signature there (Mouillet-Richard & Laurent-Puig 2020) |
| PTGER4 | Sphyris et al., Cancers (Basel), 2021 (review) | [C] C: Druggable receptor controlling YAP activity in foetal-like Sca1+ reserve-like stem cells (Sphyris 2021) |
| SNAI1 | De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Review: EMT TF driving downregulation of epithelial/stem markers in metastasis-initiating cCSCs (De Angelis 2021) |
| SNAI2 | De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Review: EMT TF in metastasis-initiating plastic cCSC state (De Angelis 2021) |
| SRC | Sphyris et al., Cancers (Basel), 2021 (review) | [B] B: SRC/YAP signalling engages foetal regenerative program (Sphyris 2021) |
| TNFAIP8 | Sphyris et al., Cancers (Basel), 2021 (review) | [B] B: Regulator of Clu+ regenerative program; KO impairs YAP nuclear recruitment post-DSS (Sphyris 2021) |
| TWIST1 | De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Review: EMT TF in metastasis-initiating plastic cCSC state (De Angelis 2021) |
| WIF1 | Ramadan et al., Trends Cancer, 2022 (review) | [B] B: Review: Wnt antagonist secreted by Apc-mutant cells to suppress normal ISCs (Ramadan 2022) |
| ZEB2 | De Angelis et al., Cancers (Basel), 2021 (review) | [C] C: Review: EMT TF in metastasis-initiating plastic cCSC state (De Angelis 2021) |

## T3 — Non-intestinal oncofetal / fetal-like references

| Oncofetal marker | 参考文献 | 注释（为何入选） |
|---|---|---|
| PLVAP | Sharma et al., Cell, 2020; Cao et al., MedComm, 2023 (review) | [A·B·C] A: Marks fetal-liver endothelial cells that reappear in tumor; cell_type endothelial, not epithelial (Sharma 2020); B: Fetal-like ecosystem re-emerges during liver regeneration; endothelial cell_type (Sharma 2020); C: Tumor-specific onco-fetal endothelial marker induced downstream of hepatocyte VEGFA; endothelial cell_type (Sharma 2020) |
| FOLR2 | Sharma et al., Cell, 2020; Cao et al., MedComm, 2023 (review) | [A·C] A: Shared marker of fetal-liver macrophages; macrophage cell_type (Sharma 2020); C: Defines fetal-like TAM1 in HCC; macrophage cell_type (Sharma 2020) |
| DCLK1 | Burdziak et al., Science, 2023; Goldenring & Mills, Gastroenterology, 2022 (review) | [B·C] B: Tuft cell marker upregulated after parietal cell loss; tuft-derived IL-25 feeds the ILC2-IL-13 axis driving SPEM (Goldenring & Mills 2022); C: Marker of the progenitor-like tuft state in premalignant pancreas; cell_type: pancreatic epithelium (Burdziak 2023) |
| IGF2BP3 | Fujiwara et al., PNAS, 2024; Viragova et al., Cell Stem Cell, 2024 (review) | [B·C] B: Fetal-like hepatocyte marker, also high in human acute liver failure (Viragova 2024); C: Termed oncofetal; correlates with high D-score; drives miR-21-5p+C isomiR and lowers Let-7 seed occupancy (Fujiwara 2024) ‖ 功能证据: siRNA KD in A549/NCI-H1703; IP with Drosha; nRIP (Fujiwara 2024) |
| SOX4 | Moorman et al., Nature, 2025; Viragova et al., Cell Stem Cell, 2024 (review) | [B·C] B: Fetal-like marker in skin wound and liver-progenitor-like cells (Viragova 2024); C: Embryonic developmental gene elevated in treatment-naive tumour cells above normal stem cells (Moorman 2025) |
| KDR | Sharma et al., Cell, 2020; Cao et al., MedComm, 2023 (review) | [C] C: Receptor defining a PLVAP+ onco-fetal endothelial subcluster; endothelial cell_type (Sharma 2020) |
| ANXA10 | Burdziak et al., Science, 2023 | [C] C: Marker of the gastric metaplastic (E6) plastic cell state in Kras-mutant pancreas; cell_type: pancreatic epithelium (Burdziak 2023) |
| CD163 | Sharma et al., Cell, 2020 | [C] C: Co-marker of fetal-like FOLR2+ TAM1; macrophage cell_type (Sharma 2020) |
| CEACAM1 | Damelin et al., Cancer Res, 2011 | [C] C: Carcinoembryonal antigen higher in 5T4hi TIC fraction (Damelin 2011) |
| CEACAM6 | Damelin et al., Cancer Res, 2011 | [C] C: Carcinoembryonal antigen higher in 5T4hi TIC fraction (Damelin 2011) |
| DLL4 | Sharma et al., Cell, 2020 | [C] C: Notch ligand on onco-fetal endothelium driving fetal-like macrophage reprogramming; endothelial cell_type (Sharma 2020) |
| ESRRB | Messier et al., J Cell Physiol, 2016 | [C] C: Example orphan nuclear receptor bivalently marked (H3K4me3+H3K27me3) in hESC, MCF7 and patient tumors; 'oncofetal epigenetic' signature; direction=bivalent, not expression (Messier 2016) |
| FABP5 | Mo et al., Hum Cell, 2023 | [C] C: Embryonic stem gene high in tumor myeloid/endothelial cells (Mo 2023) |
| HES1 | Sharma et al., Cell, 2020 | [C] C: Conserved regulon of onco-fetal macrophage reprogramming; macrophage cell_type (Sharma 2020) |
| HMGA2 | Fujiwara et al., PNAS, 2024 | [C] C: Named as an IGF2BP3 target mRNA protected from Let-7 silencing; cell_type: lung tumor epithelium (Fujiwara 2024) |
| IFITM1 | Mo et al., Hum Cell, 2023 | [C] C: Defines oncofetal CAF subtype; poor prognosis (Mo 2023) |
| IL18 | Burdziak et al., Science, 2023 | [C] C: Communication-gene marker of the gastric (E6) plastic epithelial state; cell_type: pancreatic epithelium (Burdziak 2023) |
| LHX4 | Messier et al., J Cell Physiol, 2016 | [C] C: Example stem cell transcriptional regulator bivalently marked (H3K4me3+H3K27me3) in hESC, MCF7 and patient tumors; 'oncofetal epigenetic' signature; direction=bivalent, not expression (Messier 2016) |
| LIN28B | Fujiwara et al., PNAS, 2024 | [C] C: Named as an IGF2BP3 target mRNA whose protection blocks Let-7 maturation; cell_type: lung tumor epithelium (Fujiwara 2024) |
| MAF | Sharma et al., Cell, 2020 | [C] C: Member of conserved fetal-like macrophage regulon module; macrophage cell_type (Sharma 2020) |
| MRC1 | Sharma et al., Cell, 2020 | [C] C: Protein co-marker of fetal-like TAM1 phenotype gained after Notch activation; macrophage cell_type (Sharma 2020) |
| NES | Burdziak et al., Science, 2023 | [C] C: Marks a Kras-mutant progenitor-like apex state with high epigenetic plasticity; cell_type: pancreatic epithelium (non-intestinal) (Burdziak 2023) ‖ 功能证据: scRNA/scATAC trajectory (CellRank, Palantir); Il33 knockdown shifts cells toward this state (Burdziak 2023) |
| NOTCH2 | Sharma et al., Cell, 2020 | [C] C: Receptor on fetal-like macrophages for endothelial DLL4; macrophage cell_type (Sharma 2020) |
| NR1H3 | Sharma et al., Cell, 2020 | [C] C: Transcription factor shared by fetal-liver and fetal-like tumor macrophages; macrophage cell_type (Sharma 2020) |
| NRP1 | Sharma et al., Cell, 2020 | [C] C: Co-marker of onco-fetal PLVAP+ endothelial subcluster; endothelial cell_type (Sharma 2020) |
| PODXL | Mo et al., Hum Cell, 2023 | [C] C: Defines oncofetal TEC subtype; poor prognosis (Mo 2023) |
| POU2F3 | Burdziak et al., Science, 2023 | [C] C: Marker of progenitor-like tuft cell state arising in Kras-mutant/regenerating pancreas; cell_type: pancreatic epithelium (Burdziak 2023) |
| SCGB1A1 | Damelin et al., Cancer Res, 2011 | [C] C: 20-fold higher in 5T4hi cells; lung stem cell marker (Damelin 2011) |
| SERPINH1 | Mo et al., Hum Cell, 2023 | [C] C: Embryonic stem gene heterogeneous; high in endothelial cells and fibroblasts (Mo 2023) |
| SFRP2 | Mo et al., Hum Cell, 2023 | [C] C: Defines oncofetal CAF subtype; poor prognosis (Mo 2023) |
| SMS | Mo et al., Hum Cell, 2023 | [C] C: Defines oncofetal TAM subtype with M2 polarization (Mo 2023) |
| SPIC | Sharma et al., Cell, 2020 | [C] C: Member of conserved fetal-like macrophage regulon module; macrophage cell_type (Sharma 2020) |
| SRSF7 | Mo et al., Hum Cell, 2023 | [C] C: Embryonic stem gene defining oncofetal CAF subtype; poor prognosis, ICB resistance (Mo 2023) |
| TPBG | Damelin et al., Cancer Res, 2011 | [C] C: Oncofetal antigen on tumor-initiating cells; lost on differentiation; poor prognosis in NSCLC (Damelin 2011) |
| TUBB | Mo et al., Hum Cell, 2023 | [C] C: Defines oncofetal TAM subtype (Mo 2023) |
| DDIT4 | Goldenring & Mills, Gastroenterology, 2022 (review); Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Upregulated in SPEM alongside DNA unwinding and damage-monitoring genes (Goldenring & Mills 2022) |
| IFRD1 | Goldenring & Mills, Gastroenterology, 2022 (review); Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Paligenosis gene implicated in intestinal injury response (Viragova 2024) |
| AFP | Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Liver developmental gene reactivated in fetal-like hepatocytes (non-intestinal) (Viragova 2024) |
| AQP5 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Induced SPEM marker; AQP5+ SPEM cells sit at bases of incomplete intestinal metaplasia glands (Goldenring & Mills 2022) |
| CDH17 | Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Development-associated gene in transient fetal-like hepatocyte signature (Viragova 2024) |
| FOXA2 | Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Developmental TF upregulated in fetal-like hepatocytes (Viragova 2024) |
| GATA4 | Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Embryonic TF re-expressed during regeneration (non-intestinal) (Viragova 2024) |
| GATA6 | Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: Developmental TF upregulated in fetal-like hepatocytes (Viragova 2024) |
| MUC5AC | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Marks the foveolar hyperplasia component of the injury-induced pyloric metaplasia lesion (Goldenring & Mills 2022) |
| MUC6 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Mucin defining the TFF2+/MUC6+/PGC-low SPEM lineage of pyloric metaplasia (Goldenring & Mills 2022) |
| PDX1 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Antral transcription factor whose expression confirms pyloric-type metaplasia after isthmal Ras activation (Goldenring & Mills 2022) |
| RAB7B | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: ATF3-induced lysosome and autophagy trafficking protein upregulated during paligenosis (Goldenring & Mills 2022) |
| SLC7A11 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: ROS-adaptive transporter required for chief cells to complete reprogramming into SPEM (Goldenring & Mills 2022) |
| SOX11 | Viragova et al., Cell Stem Cell, 2024 (review) | [B] B: TF driving fetal-like wound signature in skin (non-intestinal) (Viragova 2024) |
| TFF1 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Marks foveolar hyperplasia within pyloric metaplasia and normal pyloric mucosa (Goldenring & Mills 2022) |
| TFF2 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Pathognomonic SPEM marker induced when chief cells reprogram during injury-induced pyloric metaplasia (Goldenring & Mills 2022) |
| WFDC2 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] B: Induced SPEM marker, possibly unique to SPEM rather than deep antral cells (Goldenring & Mills 2022) |

## Down — 在 fetal-like / oncofetal 状态中下调的基因

| Marker (down) | 参考文献 | 注释 |
|---|---|---|
| LGR5 | Gregorieff et al., Nature, 2015; Nusse et al., Nature, 2018; Yui et al., Cell Stem Cell, 2018; Ayyaz et al., Nature, 2019; Cheung et al., Cell Stem Cell, 2020; Ganesh et al., Nat Cancer, 2020; Cañellas-Socias et al., Nature, 2022; Gil Vazquez et al., Cell Stem Cell, 2022; Heinz et al., Cancer Res, 2022; Álvarez-Varela et al., Nat Cancer, 2022; Chen et al., Cell Stem Cell, 2023; Hartl et al., Sci Adv, 2024; Moorman et al., Nature, 2025; Mzoughi et al., Nat Genet, 2025; Ogden et al., Cell Genom, 2025; van der Net et al., Cell Rep, 2025; Buissant des Amorie et al., Nature, 2026; Centonze et al., Cancer Discov, 2026; Morgan et al., Br J Cancer, 2018 (review); Beumer & Clevers, Nat Rev Mol Cell Biol, 2021 (review); De Angelis et al., Cancers (Basel), 2021 (review); Sphyris et al., Cancers (Basel), 2021 (review); Palikuqi et al., Cold Spring Harb Perspect Biol, 2022 (review); Ramadan et al., Trends Cancer, 2022 (review); Pérez-González et al., Nat Cancer, 2023 (review); Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Tape, Trends Cancer, 2024 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [A·B·C] Appears only in second-trimester mature cells, absent from first-trimester progenitors (Moorman 2025) |
| OLFM4 | Gregorieff et al., Nature, 2015; Nusse et al., Nature, 2018; Yui et al., Cell Stem Cell, 2018; Ayyaz et al., Nature, 2019; Cheung et al., Cell Stem Cell, 2020; Cañellas-Socias et al., Nature, 2022; Heinz et al., Cancer Res, 2022; Chen et al., Cell Stem Cell, 2023; Moorman et al., Nature, 2025; Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Tape, Trends Cancer, 2024 (review); Viragova et al., Cell Stem Cell, 2024 (review); Clevers, Cell, 2026 (review) | [B·C] Clu-GFP+ revival cells show no overlap with OLFM4+ CBCs (Ayyaz 2019) |
| ASCL2 | Yui et al., Cell Stem Cell, 2018; Ayyaz et al., Nature, 2019; Cañellas-Socias et al., Nature, 2022; Heinz et al., Cancer Res, 2022; Mzoughi et al., Nat Genet, 2025; van der Net et al., Cell Rep, 2025; Centonze et al., Cancer Discov, 2026; Viragova et al., Cell Stem Cell, 2024 (review) | [B·C] CBC marker of SSC2a; Clu-high revival cells rarely co-express it (Ayyaz 2019) |
| AXIN2 | Gregorieff et al., Nature, 2015; Cheung et al., Cell Stem Cell, 2020; Cañellas-Socias et al., Nature, 2022; van der Net et al., Cell Rep, 2025; Centonze et al., Cancer Discov, 2026 | [B·C] Homeostatic stem/DCS/Wnt markers lost in wound-healing state (Cheung 2020) |
| CDX2 | Moorman et al., Nature, 2025; Mzoughi et al., Nat Genet, 2025; Centonze et al., Cancer Discov, 2026; Mouillet-Richard & Laurent-Puig, Cancers (Basel), 2020 (review); Higa & Nakayama, Cancer Sci, 2024 (review) | [C] Intestinal identity TF suppressed in KRAS/MAPK-mutant HRC-high tumors (Centonze 2026) |
| SMOC2 | Cañellas-Socias et al., Nature, 2022; Heinz et al., Cancer Res, 2022; Álvarez-Varela et al., Nat Cancer, 2022; Centonze et al., Cancer Discov, 2026 | [C] ISC gene marking Emp1-TOM-low cells, i.e. low in HRCs (Cañellas-Socias 2022) |
| KRT20 | Cañellas-Socias et al., Nature, 2022; Hartl et al., Sci Adv, 2024; Moorman et al., Nature, 2025 | [B·C] Differentiation (colonocyte) marker lost in regenerative/locked-regenerative state (Hartl 2024) |
| LATS2 | Cheung et al., Cell Stem Cell, 2020; Hartl et al., Sci Adv, 2024; Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review) | [B·C] Loss of Hippo kinase induces regenerative state (Cheung 2020) |
| LRIG1 | Yui et al., Cell Stem Cell, 2018; Qin et al., Cell, 2023; Sphyris et al., Cancers (Basel), 2021 (review) | [B·C] Adult stem marker lost at protein level during repair (Yui 2018) |
| ATOH1 | Yui et al., Cell Stem Cell, 2018; Hartl et al., Sci Adv, 2024 | [B] Differentiation (secretory) marker lost in regenerative/locked-regenerative state (Hartl 2024) |
| CDX1 | Moorman et al., Nature, 2025; Centonze et al., Cancer Discov, 2026 | [C] Intestinal identity TF restored by KRAS inhibition, i.e. low in HRC state (Centonze 2026) |
| COL17A1 | Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review); Higa & Nakayama, Cancer Sci, 2024 (review) | [C] Its chemotherapy-induced loss switches on YAP activity in surviving quiescent LGR5+ tumour cells (Cañellas-Socias 2024) |
| EPHB2 | Gil Vazquez et al., Cell Stem Cell, 2022; Qin et al., Cell, 2023 | [C] Used as CBC (crypt base columnar) surface marker opposite to Ly6a; marks the non-regenerative stem pole (Gil Vazquez 2022) |
| HNF4A | Mzoughi et al., Nat Genet, 2025; Centonze et al., Cancer Discov, 2026 | [C] Intestinal identity TF suppressed in HRC/oncofetal-like KRAS-mutant tumors (Centonze 2026) |
| LATS1 | Cheung et al., Cell Stem Cell, 2020; Cañellas-Socias et al., Nat Rev Gastroenterol Hepatol, 2024 (review) | [B·C] Loss of Hippo kinase induces regenerative state (Cheung 2020) |
| LYZ | Yui et al., Cell Stem Cell, 2018; Cheung et al., Cell Stem Cell, 2020 | [B] Secretory marker absent in YAP-induced state (Cheung 2020) |
| MUC2 | Yui et al., Cell Stem Cell, 2018; Hartl et al., Sci Adv, 2024 | [B] Differentiation (goblet) marker lost in regenerative/locked-regenerative state (Hartl 2024) |
| SLC26A3 | Hartl et al., Sci Adv, 2024; Moorman et al., Nature, 2025 | [A·B] Mature enterocyte marker restricted to second-trimester colon, not fetal progenitors (Moorman 2025) |
| AQP4 | Gregorieff et al., Nature, 2015 | [B] Stem signature gene repressed by Yap in the regenerative state (Gregorieff 2015) |
| AQP8 | Hartl et al., Sci Adv, 2024 | [B] Differentiation (colonocyte) marker lost in regenerative/locked-regenerative state (Hartl 2024) |
| BHLHA15 | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] Mature chief cell transcription factor downregulated as the earliest event in SPEM reprogramming (Goldenring & Mills 2022) |
| BIRC5 | Qin et al., Cell, 2023 | [C] proCSC marker; proliferative state opposing revCSC (Qin 2023) |
| Ccl6 (mouse) | Gregorieff et al., Nature, 2015 | [B] Paneth-associated gene suppressed by Yap during regeneration (Gregorieff 2015) |
| CCN1 | Cañellas-Socias et al., Nature, 2022 | [C] Canonical YAP target not upregulated in HRCs, unlike in chemotherapy-induced fetal-like state (Cañellas-Socias 2022) |
| CCN2 | Cañellas-Socias et al., Nature, 2022 | [C] Canonical YAP target NOT upregulated in HRCs; used to show the HRC state is distinct from the YAP/fetal program (Cañellas-Socias 2022) |
| CD24 | Damelin et al., Cancer Res, 2011 | [C] Low in 5T4+ TICs; mutually exclusive with 5T4; rises on differentiation (Damelin 2011) |
| CD36 | Ganesh et al., Nat Cancer, 2020 | [C] L1CAM-high metastasis-initiating cells express low CD36; negative marker of this state (Ganesh 2020) |
| CDH1 | Ganesh et al., Nat Cancer, 2020 | [C] Loss of membrane E-cadherin contacts drives the regenerative L1CAM-high state (Ganesh 2020) |
| CHGA | Yui et al., Cell Stem Cell, 2018 | [B] Enteroendocrine marker repressed in reprogrammed fetal-like cultures (Yui 2018) |
| DACH1 | Ogden et al., Cell Genom, 2025 | [C] Stem cell marker repressed during transition into (i)REC state (Ogden 2025) |
| EPHB3 | Gregorieff et al., Nature, 2015 | [B] Wnt/stem gene repressed by Yap in the regenerative state (Gregorieff 2015) |
| GCG | Hartl et al., Sci Adv, 2024 | [B] Differentiation (enteroendocrine) marker lost in regenerative/locked-regenerative state (Hartl 2024) |
| KDM1A | Viragova et al., Cell Stem Cell, 2024 (review) | [B] Repressor of fetal/neonatal genes; decreases during radiation injury (Viragova 2024) |
| KIT | Gregorieff et al., Nature, 2015 | [B] Paneth marker suppressed by Yap during regeneration (Gregorieff 2015) |
| LCT | Yui et al., Cell Stem Cell, 2018 | [B] Enterocyte differentiation marker repressed in the fetal-like state (Yui 2018) |
| MIR148A | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] Chief-cell-enriched microRNA lost early during reprogramming (Goldenring & Mills 2022) |
| MKI67 | Rehman et al., Cell, 2021 | [C] Reduced proliferation marker in on-treatment persister tumors (Rehman 2021) |
| MYC | Rehman et al., Cell, 2021 | [C] Myc response module suppressed in diapause-like persister tumors, mirroring embryonic diapause (Rehman 2021) |
| NKD1 | Gregorieff et al., Nature, 2015 | [B] Wnt target repressed by Yap in the regenerative state (Gregorieff 2015) |
| NOTUM | Centonze et al., Cancer Discov, 2026 | [C] Stem marker gene low in HRC state, induced on KRASG12D inhibition (Centonze 2026) |
| NR5A2 | Burdziak et al., Science, 2023 | [C] Acinar identity gene defining the differentiated state that is lost with increasing plasticity; cell_type: pancreatic epithelium (Burdziak 2023) |
| PGC | Goldenring & Mills, Gastroenterology, 2022 (review) | [B] SPEM cells retain chief cell markers only at low abundance, so low PGC is part of the SPEM definition (Goldenring & Mills 2022) |
| REG4 | Cheung et al., Cell Stem Cell, 2020 | [B·C] Homeostatic stem/DCS/Wnt markers lost in wound-healing state (Cheung 2020) |
| REST | Ganesh et al., Nat Cancer, 2020 | [C] Downregulated on loss of epithelial integrity, derepressing L1CAM; entry into regenerative L1CAM-high state (Ganesh 2020) |
| RXRA | Mzoughi et al., Nat Genet, 2025 | [C] Gatekeeper lost after APC loss; its blockade induces oncofetal state (Mzoughi 2025) |
| SPDEF | Gregorieff et al., Nature, 2015 | [B] Paneth/secretory marker suppressed by Yap during regeneration (Gregorieff 2015) |
| STK3 | Cheung et al., Cell Stem Cell, 2020 | [B] Loss of Hippo kinase induces YAP-active state (Cheung 2020) |
| STK4 | Cheung et al., Cell Stem Cell, 2020 | [B] Loss of Hippo kinase induces YAP-active state (Cheung 2020) |
| TCF7 | Cheung et al., Cell Stem Cell, 2020 | [C] Only Wnt cofactor downregulated upon YAP activation (Cheung 2020) |
| TFF3 | Moorman et al., Nature, 2025 | [A] Mature goblet marker restricted to second-trimester colon, not fetal progenitors (Moorman 2025) |
| WNT3 | Gregorieff et al., Nature, 2015 | [B] Paneth marker suppressed by Yap during regeneration (Gregorieff 2015) |
